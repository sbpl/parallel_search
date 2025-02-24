#ifndef INSAT_PLANNER_MH_HPP
#define INSAT_PLANNER_MH_HPP

#include <future>
#include <utility>
#include <algorithm>
#include <iostream>
#include "planners/Planner.hpp"
#include <common/insat/InsatState.hpp>
#include <common/insat/InsatEdge.hpp>
#include <common/State.hpp>

namespace ps
{

    class InsatPlannerMH : virtual public Planner
    {
    public:
        int hidx_;
        double w2_;
        // Typedefs
        typedef std::unordered_map<size_t, InsatStatePtrType> InsatStatePtrMapType;
        // typedef smpl::intrusive_heap<InsatState, IsLesserState> InsatStateQueueMinType;
        typedef smpl::intrusive_heap<InsatState::HeapData, IsLesserHeapData> InsatStateQueueMinType;

        InsatPlannerMH(ParamsType planner_params);;

        ~InsatPlannerMH() {};

        void SetStartState(const StateVarsType& state_vars);

        bool Plan();

        TrajType getSolutionTraj();
        double computeAncHeuristic(const StatePtrType& state_ptr);

    protected:
        void initialize();

        std::vector<InsatStatePtrType> getStateAncestors(const InsatStatePtrType state_ptr, bool reverse=false) const;

        void expandState(InsatStatePtrType state_ptr);

        void updateState(InsatStatePtrType& state_ptr,
                         std::vector<InsatStatePtrType>& ancestors,
                         InsatActionPtrType& action_ptr,
                         ActionSuccessor& action_successor);

        void constructInsatActions();
        void RemoveFromAllInadOpenLists(InsatState* state);
        void RemoveFromAncOpenList(InsatState* state);

        InsatStatePtrType constructInsatState(const StateVarsType& state);

        void cleanUp();

        void resetStates();

        void constructPlan(InsatStatePtrType& insat_state_ptr);

        void exit();


        std::vector<std::shared_ptr<InsatAction>> insat_actions_ptrs_;
        InsatStatePtrType start_state_ptr_;
        InsatStatePtrType goal_state_ptr_;
        InsatStateQueueMinType insat_state_open_list_;
        // vector for multiple heuristics
        std::vector<InsatStateQueueMinType> insat_open_lists_;
        // InsatStateQueueMinType anchor_list_;
        InsatStatePtrMapType insat_state_map_;
        TrajType soln_traj_;

    };

}

#endif
