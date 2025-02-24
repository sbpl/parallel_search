/*
 * Copyright (c) 2023, Ramkumar Natarajan
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Carnegie Mellon University nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
/*!
 * \file   InsatPlanner.cpp
 * \author Ramkumar Natarajan (rnataraj@cs.cmu.edu)
 * \date   2/26/23
 */

#include <planners/insat/InsatPlannerMH.hpp>
#include <vector>                                          // For std::vector
#include <Eigen/Dense>                                     // For MatrixX
#include <drake/common/trajectories/bspline_trajectory.h>  // For BsplineTrajectory
#include <drake/math/bspline_basis.h>                     // For BsplineBasis

namespace ps
{


    InsatPlannerMH::InsatPlannerMH(ParamsType planner_params) :
            Planner(planner_params)
    {
        if (planner_params.find("adaptive_opt") == planner_params.end())
        {
            planner_params["adaptive_opt"] = false;
        }
    }

    void InsatPlannerMH::SetStartState(const StateVarsType &state_vars) {
        start_state_ptr_ = constructInsatState(state_vars);
        
    }

    // NEW FUNCTIONS

    InsatState* GetInsatStateFromHeapData(const State::HeapData* d)
    {
        return (InsatState*)((char*)d - d->off);
    }

    void InsatPlannerMH::RemoveFromAllInadOpenLists(InsatState* state)
    {
        for (size_t hidx = 1; hidx < num_heuristics_ + 1; ++hidx) {
            if (insat_open_lists_[hidx].contains(&state->open_data[hidx])) {
                insat_open_lists_[hidx].erase(&state->open_data[hidx]);
            }
        }
    }

    void InsatPlannerMH::RemoveFromAncOpenList(InsatState* state)
    {
        if (insat_open_lists_[0].contains(&state->open_data[0])) {
            insat_open_lists_[0].erase(&state->open_data[0]);
        }
    }

    bool InsatPlannerMH::Plan() {
        initialize();
        startTimer();
        while (!insat_open_lists_[0].empty() && !checkTimeout())
        {
            // instead of popping min, pop based on heuristic value  
            hidx_ %= num_heuristics_;
            hidx_ ++;

            InsatState* state_ptr;

            auto anc_min = insat_open_lists_[0].min();
            if (!insat_open_lists_[hidx_].empty()){
                auto open_min = insat_open_lists_[hidx_].min();

                if (open_min->f <= w2_*anc_min->f){
                    // pick state from open list
                    state_ptr = GetInsatStateFromHeapData(insat_open_lists_[hidx_].min());
                    state_ptr->is_visited_mh = true;
                    RemoveFromAllInadOpenLists(state_ptr);
                }
                else{
                    // pick state from anchor list
                    state_ptr = GetInsatStateFromHeapData(insat_open_lists_[0].min());
                    state_ptr->is_visited_anc = true;
                    RemoveFromAncOpenList(state_ptr);
                }
            }
            else{
                // pick state from anchor list
                state_ptr = GetInsatStateFromHeapData(insat_open_lists_[0].min());
                state_ptr->is_visited_anc = true;
                RemoveFromAncOpenList(state_ptr);
            }
            // auto state_ptr = insat_open_lists[hidx].min();
            // insat_state_open_list_.pop();
            // save state ptr
            state_ptrs_all_.push_back(state_ptr);

            // Return solution if goal state is expanded
            if (isGoalState(state_ptr))
            {
                auto t_end = std::chrono::steady_clock::now();
                double t_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_end-t_start_).count();
                goal_state_ptr_ = state_ptr;

                // Reconstruct and return path
                constructPlan(state_ptr);
                planner_stats_.total_time_ = 1e-9*t_elapsed;
                exit();
                return true;
            }

            // set closed in all lists 
            expandState(state_ptr);

        }

        auto t_end = std::chrono::steady_clock::now();
        double t_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_end-t_start_).count();
        planner_stats_.total_time_ = 1e-9*t_elapsed;
        return false;
    }

    TrajType InsatPlannerMH::getSolutionTraj() {
        return soln_traj_;
    }

    double InsatPlannerMH::computeAncHeuristic(const StatePtrType& state_ptr){
        double anc_heur = 1000000;
        for (auto goal : goals_list_){
            StatePtrType goal_state_ptr = new InsatState(goal);
            auto heur = computeHeuristic(state_ptr, goal_state_ptr);
            anc_heur = std::min(heur, anc_heur);
        }
        return anc_heur;
    }

    void InsatPlannerMH::initialize() {

        plan_.clear();
        hidx_ = 0;
        w2_ = 1.3;
        insat_open_lists_.resize(num_heuristics_+1);

        // Initialize planner stats
        planner_stats_ = PlannerStats();

        // Initialize start state
        start_state_ptr_->SetGValue(0);
        // start_state_ptr_->SetHValue(computeHeuristic(start_state_ptr_));

        double h_val = computeAncHeuristic(start_state_ptr_); // Compute heuristic for each index
        start_state_ptr_->open_data[0].h = h_val;
        insat_open_lists_[0].push(&start_state_ptr_->open_data[0]);

        // std::cout << "GOALS LIST IS" << goals_list_.size() << std::endl;

        for (int i = 0; i < num_heuristics_; ++i) { //check this + 1 logic
            // need to convert this goal into a state pointer type to feed into compute heuristic
            StatePtrType goal_state_ptr = new InsatState(goals_list_[i]);

            double h_val = computeHeuristic(start_state_ptr_, goal_state_ptr);

            // std::cout << "After computing heuristic:" << std::endl;
            // start_state_ptr_->Print("DEBUG: ");

            // h_val = computeHeuristic(start_state_ptr_, goal_state_ptr); // Compute heuristic for each index
            start_state_ptr_->open_data[i+1].h = h_val;
 
            double f_val = start_state_ptr_->GetGValue() + heuristic_w_*h_val;
            start_state_ptr_->open_data[i+1].f = f_val;
            // start_state_ptr_->SetFValue(f_val);  heuristic_weights_[i] * (weights?)

            insat_open_lists_[i+1].push(&start_state_ptr_->open_data[i]);
        }
        // for (size_t hidx = 0; hidx < NumHeuristics(search); ++hidx) {
        //     search.open_lists_[hidx].push(&search.start_state_->open_data[hidx]);
        //     SMPL_DEBUG("Inserted start state %d into search %zu with f = %d",
        //             start_id,
        //             hidx,
        //             ComputeFVal(search, search.start_state_->open_data[hidx]));
        // }

        // anchor_list_ = insat_open_lists_[0];

        // Reset goal state
        goal_state_ptr_ = NULL;

        // Reset state
        planner_stats_ = PlannerStats();

        // Reset h_min
        h_val_min_ = DINF;

        planner_stats_.num_jobs_per_thread_.resize(1, 0);
        // Initialize open list
        // start_state_ptr_->SetFValue(start_state_ptr_->GetGValue() + heuristic_w_*start_state_ptr_->GetHValue());
        // insat_state_open_list_.push(start_state_ptr_);

        constructInsatActions();
    }

    std::vector<InsatStatePtrType>
    InsatPlannerMH::getStateAncestors(const InsatStatePtrType state_ptr, bool reverse) const {
        // Get ancestors
        std::vector<InsatStatePtrType> ancestors;
        ancestors.push_back(state_ptr);
        auto bp = state_ptr->GetIncomingInsatEdgePtr();
        while (bp)
        {
            ancestors.push_back(bp->lowD_parent_state_ptr_);
            bp = bp->lowD_parent_state_ptr_->GetIncomingInsatEdgePtr();
        }
        if (reverse)
        {
            std::reverse(ancestors.begin(), ancestors.end());
        }
        return ancestors;
    }

    void InsatPlannerMH::expandState(InsatStatePtrType state_ptr) {

        if (VERBOSE) state_ptr->Print("Expanding");
        planner_stats_.num_state_expansions_++;

        // state_ptr->SetVisited();

        auto ancestors = getStateAncestors(state_ptr);

        for (auto& action_ptr: insat_actions_ptrs_)
        {
            if (action_ptr->CheckPreconditions(state_ptr->GetStateVars()))
            {
                // Evaluate the edge
                auto action_successor = action_ptr->GetSuccessor(state_ptr->GetStateVars());

                updateState(state_ptr, ancestors, action_ptr, action_successor);
            }
        }
    }

    void InsatPlannerMH::updateState(InsatStatePtrType &state_ptr, std::vector<InsatStatePtrType> &ancestors,
                                   InsatActionPtrType &action_ptr, ActionSuccessor &action_successor) {
        planner_stats_.num_evaluated_edges_++;

        if (action_successor.success_)
        {
            auto successor_state_ptr = constructInsatState(action_successor.successor_state_vars_costs_.back().first);

            if (!successor_state_ptr->IsVisitedAnc())
            {
                InsatStatePtrType best_anc;
                TrajType traj;
                double cost = 0;
                double inc_cost = 0;

                if (planner_params_["smart_opt"] == true)
                {
                    std::vector<StateVarsType> anc_states;
                    for (auto& anc: ancestors)
                    {
                        anc_states.emplace_back(anc->GetStateVars());
                    }

                    if (state_ptr->GetIncomingInsatEdgePtr()) /// When anc is not start
                    {
                        traj = action_ptr->optimize(state_ptr->GetIncomingInsatEdgePtr()->GetTraj(),
                                                    anc_states,
                                                    successor_state_ptr->GetStateVars());
                        inc_cost = action_ptr->getCost(traj) - action_ptr->getCost(state_ptr->GetIncomingInsatEdgePtr()->GetTraj());
                    }
                    else
                    {
                        traj = action_ptr->optimize(TrajType(),
                                                    anc_states,
                                                    successor_state_ptr->GetStateVars());
                        inc_cost = action_ptr->getCost(traj);
                    }
                    if (traj.isValid())
                    {
                        best_anc = start_state_ptr_;
                    }
                }
                else
                {
                    for (auto& anc: ancestors)
                    {
                        if (planner_params_["adaptive_opt"] == true)
                        {
                            if (anc->GetIncomingInsatEdgePtr()) /// When anc is not start
                            {
                                traj = action_ptr->optimize(anc->GetIncomingInsatEdgePtr()->GetTraj(),
                                                            anc->GetStateVars(),
                                                            successor_state_ptr->GetStateVars());
                                inc_cost = action_ptr->getCost(traj) - action_ptr->getCost(
                                    anc->GetIncomingInsatEdgePtr()->GetTraj());
                            }
                            else
                            {
                                traj = action_ptr->optimize(TrajType(),
                                                            anc->GetStateVars(),
                                                            successor_state_ptr->GetStateVars());
                                inc_cost = action_ptr->getCost(traj);
                            }
                        }
                        else
                        {
                            TrajType inc_traj = action_ptr->optimize(anc->GetStateVars(), successor_state_ptr->GetStateVars());
                            if (inc_traj.size() > 0)
                            {
                                inc_cost = action_ptr->getCost(inc_traj);
                                if (anc->GetIncomingInsatEdgePtr()) /// When anc is not start
                                {
                                    traj = action_ptr->warmOptimize(anc->GetIncomingInsatEdgePtr()->GetTraj(), inc_traj);
                                }
                                else
                                {
                                    traj = action_ptr->warmOptimize(inc_traj);
                                }
                            }
                            else
                            {
                                continue;
                            }

                        }
                        if (traj.isValid())
                        {
                            best_anc = anc;
                        }
                    }
                }

                if (traj.disc_traj_.cols()<=2)
                {
                    return;
                }

                cost = action_ptr->getCost(traj);
                double new_g_val = cost;

                if (successor_state_ptr->GetGValue() - new_g_val > 0.01)
                {

                    // double h_val = successor_state_ptr->GetHValue(); //replace with getting ith h value
                    // if (h_val == -1)
                    // {
                    //     h_val = computeHeuristic(successor_state_ptr);  // change to computing distance to i'th goal
                    //     successor_state_ptr->SetHValue(h_val); // i'th
                    // }

                    // if (h_val != DINF)
                    // {
                    //     h_val_min_ = h_val < h_val_min_ ? h_val : h_val_min_;   // where is this coming from?
                    //     successor_state_ptr->SetGValue(new_g_val); //
                    //     successor_state_ptr->SetFValue(new_g_val + heuristic_w_*h_val); //

                    //     auto edge_ptr = new Edge(state_ptr, action_ptr, successor_state_ptr);
                    //     edge_ptr->SetCost(inc_cost);
                    //     successor_state_ptr->SetIncomingEdgePtr(edge_ptr);

                    //     auto insat_edge_ptr = new InsatEdge(state_ptr, action_ptr, best_anc, successor_state_ptr);
                    //     insat_edge_ptr->SetTraj(traj);
                    //     insat_edge_ptr->SetTrajCost(cost);
                    //     insat_edge_ptr->SetCost(cost);
                    //     if (isGoalState(successor_state_ptr))
                    //     {
                    //         insat_edge_ptr->SetTrajCost(0);
                    //         insat_edge_ptr->SetCost(0);
                    //         successor_state_ptr->SetFValue(0.0);
                    //     }
                    //     edge_map_.insert(std::make_pair(getEdgeKey(insat_edge_ptr), insat_edge_ptr));
                    //     successor_state_ptr->SetIncomingInsatEdgePtr(insat_edge_ptr); //

                    //     if (insat_state_open_list_.contains(successor_state_ptr))  //how to deal with this in a multi-heuristic case? 
                    //     {
                    //         insat_state_open_list_.decrease(successor_state_ptr);
                    //     }
                    //     else
                    //     {
                    //         insat_state_open_list_.push(successor_state_ptr);
                    //     }
                    // }

                    auto edge_ptr = new Edge(state_ptr, action_ptr, successor_state_ptr);
                    edge_ptr->SetCost(inc_cost);
                    successor_state_ptr->SetIncomingEdgePtr(edge_ptr);

                    auto insat_edge_ptr = new InsatEdge(state_ptr, action_ptr, best_anc, successor_state_ptr);
                    insat_edge_ptr->SetTraj(traj);
                    insat_edge_ptr->SetTrajCost(cost);
                    insat_edge_ptr->SetCost(cost);
                    if (isGoalState(successor_state_ptr))
                    {
                        insat_edge_ptr->SetTrajCost(0);
                        insat_edge_ptr->SetCost(0);
                        successor_state_ptr->SetFValue(0.0);
                    }
                    edge_map_.insert(std::make_pair(getEdgeKey(insat_edge_ptr), insat_edge_ptr));
                    successor_state_ptr->SetIncomingInsatEdgePtr(insat_edge_ptr); //

                    // update f anc
                    double h_val = computeAncHeuristic(successor_state_ptr); // Compute heuristic for each index
                    successor_state_ptr->open_data[0].h = h_val;
                    successor_state_ptr->SetGValue(new_g_val); //
                    successor_state_ptr->open_data[0].f = new_g_val + heuristic_w_*h_val; //
                    // push to anc?
                    // insat_open_lists_[0].push(&successor_state_ptr->open_data[0]);

                    if (!successor_state_ptr->IsVisitedInad()){
                        for (int i = 0; i < num_heuristics_; ++i) {
                            // need to convert this goal into a state pointer type to feed into compute heuristic
                            StatePtrType goal_state_ptr = new InsatState(goals_list_[i]);
                            h_val = computeHeuristic(successor_state_ptr, goal_state_ptr); // Compute heuristic for each index
                            successor_state_ptr->open_data[i+1].h = h_val;
                            successor_state_ptr->open_data[i+1].f = new_g_val + heuristic_w_*h_val; //
                            // start_state_ptr_->SetFValue(f_val);  heuristic_weights_[i] * (weights?)

                            // insat_open_lists_[i+1].push(&successor_state_ptr->open_data[i+1]);
                        }
                    }

                    

                //     if (insat_state_open_list_.contains(successor_state_ptr->open_data))  //how to deal with this in a multi-heuristic case? 
                //     {
                //         insat_state_open_list_.decrease(successor_state_ptr);
                //     }
                //     else
                //     {
                //         insat_state_open_list_.push(successor_state_ptr);
                //     }
                // }
                    for (size_t i = 0; i < num_heuristics_+1; ++i) {
                        if (insat_open_lists_[i].contains(&successor_state_ptr->open_data[i])) {
                            // Decrease key for the specific heuristic's queue.
                            insat_open_lists_[i].decrease(&successor_state_ptr->open_data[i]);
                        } else {
                            // Push the state into the specific heuristic's queue.
                            insat_open_lists_[i].push(&successor_state_ptr->open_data[i]);
                        }
                    }
                }
            }
        }
    }

    void InsatPlannerMH::constructInsatActions() {
        for (auto& action_ptr : actions_ptrs_)
        {
            insat_actions_ptrs_.emplace_back(std::dynamic_pointer_cast<InsatAction>(action_ptr));
        }
    }

    InsatStatePtrType InsatPlannerMH::constructInsatState(const StateVarsType &state) {
        size_t key = state_key_generator_(state);
        auto it = insat_state_map_.find(key);
        InsatStatePtrType insat_state_ptr;

        // Check if state exists in the search state map
        if (it == insat_state_map_.end())
        {
            insat_state_ptr = new InsatState(state);
            insat_state_ptr->open_data.resize(num_heuristics_ + 1);
            for (size_t i = 0; i < num_heuristics_+1; ++i) {
                // insat_state_ptr->open_data[i].off = ((char*)&insat_state_ptr->open_data[0] - (char*)insat_state_ptr);
                // insat_state_ptr->open_data[i].h = -1;
                insat_state_ptr->open_data.resize(num_heuristics_ + 1);
                for (size_t i = 0; i < num_heuristics_ + 1; ++i) {
                    State::HeapData* d = &insat_state_ptr->open_data[i];  // Get HeapData pointer
                    d->off = (char*)d - (char*)insat_state_ptr;    // Store the correct offset
                    d->h = -1;
                }

            }
            insat_state_map_.insert(std::pair<size_t, InsatStatePtrType>(key, insat_state_ptr));
        }
        else
        {
            insat_state_ptr = it->second;
        }

        return insat_state_ptr;
    }

    void InsatPlannerMH::cleanUp() {
        for (auto& state_it : insat_state_map_)
        {
            if (state_it.second)
            {
                delete state_it.second;
                state_it.second = NULL;
            }
        }
        insat_state_map_.clear();

        for (auto& edge_it : edge_map_)
        {
            if (edge_it.second)
            {
                delete edge_it.second;
                edge_it.second = NULL;
            }
        }
        edge_map_.clear();

        State::ResetStateIDCounter();
        Edge::ResetStateIDCounter();
    }

    void InsatPlannerMH::resetStates() {
        for (auto it = insat_state_map_.begin(); it != insat_state_map_.end(); ++it)
        {
            it->second->ResetGValue();
            // it->second->ResetFValue();
            // it->second->ResetVValue();
            it->second->ResetIncomingInsatEdgePtr();
            // it->second->UnsetVisited();
            it->second->UnsetBeingExpanded();
            it->second->num_successors_ = 0;
            it->second->num_expanded_successors_ = 0;
        }
    }

    void InsatPlannerMH::constructPlan(InsatStatePtrType &insat_state_ptr) {
      StatePtrType state_ptr = insat_state_ptr;
      Planner::constructPlan(state_ptr);

      if (insat_state_ptr->GetIncomingInsatEdgePtr())
        {
//                planner_stats_.path_cost_ = insat_state_ptr->GetIncomingInsatEdgePtr()->GetTrajCost();
            planner_stats_.path_cost_ =
                    insat_actions_ptrs_[0]->getCost(insat_state_ptr->GetIncomingInsatEdgePtr()->GetTraj());
            soln_traj_ = insat_state_ptr->GetIncomingInsatEdgePtr()->GetTraj();
        }
    }

    void InsatPlannerMH::exit() {
        // Clear open list
        // while (!insat_state_open_list_.empty())
        // {
        //     insat_state_open_list_.pop();
        // }

        cleanUp();
    }
}
