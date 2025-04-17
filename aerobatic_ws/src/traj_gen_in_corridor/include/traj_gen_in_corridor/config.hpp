#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <vector>

#include <ros/ros.h>

struct Config
{
    std::string targetTopic, triggerTopic;
    std::string odomTopic;
    std::string trajTopic, visTrajTopic;
    double maxOmgRate;
    double maxVelRate;
    double maxTau, minTau;
    double maxRange;
    double weightT;
    std::vector<double> chiVec;
    double smoothingEps;
    double relCostTol;
    double quadratureResolution;
    double flipResolution;
    std::vector<double> set_att;
    std::vector<double> set_time;
    double deltaT;
    bool isDebug;
    bool isOpt_settime;

    // Load all parameters specified by ROS script
    inline void load(const ros::NodeHandle &nh_priv)
    {
        nh_priv.getParam("TargetTopic", targetTopic);
        nh_priv.getParam("TriggerTopic", triggerTopic);
        nh_priv.getParam("OdomTopic", odomTopic);
        nh_priv.getParam("TrajTopic", trajTopic);
        nh_priv.getParam("VisTrajTopic", visTrajTopic);
        nh_priv.getParam("MaxOmgRate", maxOmgRate);
        nh_priv.getParam("MaxVelRate", maxVelRate);
        nh_priv.getParam("MaxTau", maxTau);
        nh_priv.getParam("MinTau", minTau);
        nh_priv.getParam("MaxRange", maxRange);
        nh_priv.getParam("WeightT", weightT);
        nh_priv.getParam("ChiVec", chiVec);
        nh_priv.getParam("SmoothingEps", smoothingEps);
        nh_priv.getParam("RelCostTol", relCostTol);
        nh_priv.getParam("QuadratureResolution", quadratureResolution);
        nh_priv.getParam("FlipResolution", flipResolution);
        nh_priv.getParam("set_att", set_att);
        nh_priv.getParam("set_time", set_time);
        nh_priv.getParam("deltaT", deltaT);
        nh_priv.getParam("debug", isDebug);
        nh_priv.getParam("opt_settime", isOpt_settime);
    }
};

#endif