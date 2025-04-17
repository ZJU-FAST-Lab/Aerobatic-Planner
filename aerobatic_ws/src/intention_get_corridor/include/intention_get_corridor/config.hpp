#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <vector>

#include <ros/ros.h>

struct Config
{
    std::string infmapTopic;            
    std::string odomTopic;              
    std::string targetTopic;           
    std::string pointCloudPath;       
    std::string corridorPath;           
    double dilateRadius;               
    double gridResolution;              
    std::vector<double> r3Bound;       
    double localBoxHalfWidth;           
    std::vector<double> expectedHeight; 
    int outlierThreshold;            
    bool useLoadPCDFile;                
    double astar_weight;               
    int cnt_pos;                        
    std::vector<double> set_att;        
    std::vector<double> set_pos;        

    inline void
    load(const ros::NodeHandle &nh_priv)
    {
        nh_priv.getParam("InfMapTopic", infmapTopic);
        nh_priv.getParam("OdomTopic", odomTopic);
        nh_priv.getParam("TargetTopic", targetTopic);
        nh_priv.getParam("PointCloudPath", pointCloudPath);
        nh_priv.getParam("CorridorPath", corridorPath);
        nh_priv.getParam("DilateRadius", dilateRadius);
        nh_priv.getParam("GridResolution", gridResolution);
        nh_priv.getParam("R3Bound", r3Bound);
        nh_priv.getParam("LocalBoxHalfWidth", localBoxHalfWidth);
        nh_priv.getParam("ExpectedHeight", expectedHeight);
        nh_priv.getParam("OutlierThreshold", outlierThreshold);
        nh_priv.getParam("PointCloudUsePCD", useLoadPCDFile);
        nh_priv.getParam("Astar_weight", astar_weight);
        nh_priv.getParam("SetPos", set_pos);
        nh_priv.getParam("SetAtt", set_att);
        cnt_pos = int(set_pos.size() / 3);
        if (int(set_pos.size()) != cnt_pos * 3 || int(set_att.size()) != cnt_pos * 4)
        {
            ROS_ERROR("set_pos count is wrong!");
            cnt_pos = 0;
        }
    }
};

#endif