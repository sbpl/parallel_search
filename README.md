# PINSAT: Parallelized Interleaving of Search and Trajectory Optimization

This software equips INSAT [[1]](https://arxiv.org/abs/2101.12548), [[2]](https://arxiv.org/abs/2210.08627) with parallel processing capability, namely [_PINSAT_](https://arxiv.org/abs/2401.08948), using the ideas developed in [ePA\*SE](https://arxiv.org/abs/2203.01369). 

This particular fork catkinizes the workspace to be able to build with ROS and implements [wrapper functions](https://github.com/sbpl/parallel_search/blob/catkin/examples/manipulation/MoveitInterface.hpp) to be used by the [Shield Planner](https://github.com/sbpl/shield_planner/)

### Dependencies
The planning code only uses C++ standard libraries and has no external dependencies. However, the robot_nav_2d example has the following dependencies:
1. [Boost C++](https://www.boost.org/)
2. OpenCV C++ for visualization
3. [sbpl package](https://github.com/sbpl/sbpl)


For the readme of pure search algorithms in this repository, please visit [here](https://github.com/shohinm/parallel_search).

If you use PINSAT, please cite
```
@article{natarajan2024pinsat,
  title={PINSAT: Parallelized Interleaving of Graph Search and Trajectory Optimization for Kinodynamic Motion Planning},
  author={Natarajan, Ramkumar and Mukherjee, Shohin and Choset, Howie and Likhachev, Maxim},
  journal={IEEE/RSJ International Conference on Intelligent Robots and Systems},
  year={2024}
}
```


