# SFUISE (Spline Fusion-Based Ultra-Wideband-Inertial State Estimation)
## -- Continuous-Time Ultra-Wideband-Inertial Fusion
This is the code repository of SFUISE, the first online continuous-time UWB-inertial state estimation system in ROS. Quaternion-based cubic B-splines are exploited for representing states continuously over time with efficient solutions to time derivatives and spatial differentiations in closed form. Based thereon, the functional core of the system, a novel sliding-window spline fitting scheme, is established with a customized implementation of LM method. 

The system supports UWB-inertial fusion for both ToA and TDoA principles of ultra-wideband ranging with minimized external dependencies. See below for usage of SFUISE using public and our own recorded data sets.

The work has been submitted to IEEE Robotics and Automation Letters (RA-L) and is currently under review. 

![eval_util](https://github.com/KIT-ISAS/SFUISE/blob/main/doc/util_sequences.gif)
## BibTex Citation
```
@article{arXiv23_Li,
 title = {{Continuous-Time Ultra-Wideband-Inertial Fusion}},
 author = {Kailai Li and Ziyu Cao and Uwe D. Hanebeck},
 journal = {},
 month = {January},
 url = {},
 year = {2023}
}
```

## Dependency
System dependencies (tested on Ubuntu 20.04)
* [ROS](http://wiki.ros.org/noetic/Installation) (tested with Noetic)
* [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) (tested with Eigen 3.3.7)
* [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html)
## Compilation
Compile with [catkin_tools](https://catkin-tools.readthedocs.io/en/latest/index.html):
```
cd ~/catkin_ws/src
git clone https://github.com/KIT-ISAS/SFUISE
cd ..
catkin build sfuise
```
## Usage
Run following commands in terminal
* Example for running sfuise on [UTIL](https://utiasdsl.github.io/util-uwb-dataset/) (TDoA-inertial):
```
# Change anchor_path in config_test_util.yaml
roslaunch sfuise sfuise_test_util.launch
rosbag play const1-trial1-tdoa2.bag
```
* Example for running sfuise on [ISAS-Walk](https://github.com/KIT-ISAS/SFUISE/tree/main/dataset) (ToA-inertial, own data set):
```
roslaunch sfuise sfuise_test_isas-walk1.launch
rosbag play ISAS-Walk1.bag
```
## Contributors
Kailai Li (Email: kailai.li@liu.se)

Ziyu Cao (Email: ziyu.cao@kit.edu)
## Credits
We hereby recommend reading [lie-spline-experiments](https://gitlab.com/tum-vision/lie-spline-experiments) for reference. The IMU integration was derived from [VINS-Fusion](https://github.com/HKUST-Aerial-Robotics/VINS-Fusion).
## License
The source code is released under [GPLv3](https://www.gnu.org/licenses/) license.

We are constantly working on improving our code. For any technical issues, please contact 
Kailai Li (kailai.li@liu.se).
