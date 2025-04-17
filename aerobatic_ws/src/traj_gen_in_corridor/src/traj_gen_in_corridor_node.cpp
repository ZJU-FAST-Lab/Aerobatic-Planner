#include "traj_gen_in_corridor/traj_gen_in_corridor.hpp"

int main(int argc, char **argv)
{
    ros::init(argc, argv, "traj_gen_in_corridor_node");
    ros::NodeHandle nh_;

    // Load the parameters
    Config config;
    config.load(ros::NodeHandle("~"));

    // Create a planner the run the visualization at a certain rate
    GlobalPlanner traj_gen_in_corridor(config, nh_);

    // traj_gen_in_corridor.loadCorridor();

    ros::spin();

    return 0;
}
