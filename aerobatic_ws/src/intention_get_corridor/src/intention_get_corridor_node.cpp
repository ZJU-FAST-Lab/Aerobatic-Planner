#include "intention_get_corridor/intention_get_corridor.hpp"

int main(int argc, char **argv)
{
    ros::init(argc, argv, "intention_get_corridor_node");
    ros::NodeHandle nh_;

    Config config;
    config.load(ros::NodeHandle("~"));

    GlobalPlanner intention_get_corridor(config, nh_);
    intention_get_corridor.initializeMap();

    ros::Rate lr(1000);
    while (ros::ok())
    {
        ros::spinOnce();
        lr.sleep();
    }

    return 0;
}
