# Aerobatic-Planner
An autonomous aerobatic system that is capable of complex flight maneuvers, typically achievable only by professional pilots.

**Paper**: [Unlocking Aerobatic Potential of Quadcopters: Autonomous Freestyle Flight Generation and Execution](https://www.science.org/doi/10.1126/scirobotics.adp9905). Science Robotics.

**Authors**: Mingyang Wang, Qianhao Wang, Ze Wang, Yuman Gao, Jingping Wang, Can Cui, Yuan Li, Ziming Ding, Kaiwei Wang, Chao Xu and [Fei Gao](http://zju-fast.com/research-group/fei-gao/) form the ZJU Fast Lab.

**Video Links**: [bilibili](https://www.bilibili.com/video/BV1WXouYmEaW/)

This work was born out of [MINCO](https://github.com/ZJU-FAST-Lab/GCOPTER). If you find it interesting, please give both repos stars generously. Thanks.

## Setting up the Development Environment

The provided code is intended to run on an Ubuntu 20.04 system, with the additional requirements of ROS (Robot Operating System) and Unity Editor (optional but recommended).

### Installing Ubuntu 20.04

Please install Ubuntu 20.04 on your computer. For installation instructions, please refer to the [official Ubuntu documentation](http://releases.ubuntu.com/20.04/).

### Installing ROS Noetic

To use our code, you will need to install ROS Noetic. For installation instructions, please refer to the [ROS noetic installation guide](http://wiki.ros.org/noetic/Installation/Ubuntu).

Make sure to install the `ros-noetic-desktop-full` package to have all the necessary components.

### Installing Unity Editor

To set up Unity for running the provided code, you will need to install Unity Hub. Unity Hub will allow you to manage various versions of Unity Editor and select the specific version required for testing our code.

Please refer to the [official Unity documentation](https://docs.unity3d.com/hub/manual/InstallHub.html#install-hub-linux) for instructions on how to install Unity Hub and Unity Editor. When installing Unity Hub, it will automatically prompt you to install a Unity Editor. Simply choose the version 2022.3.21f1 to align with the requirements of our code.

Once Ubuntu 20.04, ROS Noetic, Unity Hub, and Unity Editor are installed, the required environment will be in place to run the provided code.

## Run the Code

Please watch the videos in `WatchMeToRunTheCode/` folder to start the code.

## Trouble Shooting

If you encounter package dependency issues during the `catkin_make` process, please execute the following command:

```bash
sudo apt-get install ros-noetic-<missing_package>
```

This will help resolve any missing package dependencies.

Another potential issue could arise if you are unable to run the `unity_sim.launch` file. This may be due to the Python script not having executable permissions. To address this, please run the following command:

```bash
sudo chmod +x <path-to-AerobaticFlight>/aerobatic_ws/src/uav_simulator/Unity_utils/ROS-TCP-Endpoint/src/ros_tcp_endpoint/default_server_endpoint.py
```

Here, `path-to-AerobaticFlight` should be replaced with the absolute path of the repository package provided.

## Licence

The source code is released under [GPLv3](https://www.gnu.org/licenses/) license.

## Maintenance

We are still working on extending the proposed system and improving code reliability.

For any technical issues, please contact Mingyang Wang (wmyang@zju.edu.cn) or Fei Gao (fgaoaa@zju.edu.cn).

For commercial inquiries, please contact Fei Gao (fgaoaa@zju.edu.cn).