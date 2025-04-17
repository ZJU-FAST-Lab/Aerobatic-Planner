#ifndef GLBMAP_HPP
#define GLBMAP_HPP

#include "intention_get_corridor/config.hpp"
#include "intention_get_corridor/dilater.hpp"
#include "intention_get_corridor/solver/quickhull.hpp"
#include "intention_get_corridor/solver/firi.hpp"
#include "intention_get_corridor/solver/geoutils.hpp"

#include <cstdint>
#include <memory>
#include <vector>
#include <random>

#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>

#include <Eigen/Eigen>

// An implementation of occupancy gridmap with
// optimal dilation complexity
class OccupancyGridMap
{
private:
    // discrete size in xyz three directions
    const int sizeX, sizeY, sizeZ;
    // total number of voxels
    const int sizeVol;
    // origin of the map
    const double oX, oY, oZ;
    // step of y and z direction in an array
    const int stepY, stepZ;
    // scale is the size of each voxel, invSacle the inverse
    const double scale, invScale;
    // true center of the voxel (0,0,0)
    const double oXc, oYc, oZc;
    // data block that stores the voxel states
    uint8_t *occPtr;
    // bound for each direction considering the step
    const int boundX, boundY, boundZ;
    // this is just for efficiency
    const double stepScaleX, stepScaleY, stepScaleZ;
    // this stores the surface of obstacles after dilation
    // obtained for free
    std::vector<Eigen::Vector3i> dilatedSurface;

public:
    OccupancyGridMap(const Eigen::Vector3i &xyz,
                     const Eigen::Vector3d &offset,
                     const double &voxelScale)
        : sizeX(xyz(0)), sizeY(xyz(1)), sizeZ(xyz(2)),
          sizeVol(sizeX * sizeY * sizeZ),
          oX(offset(0)), oY(offset(1)), oZ(offset(2)),
          stepY(sizeX), stepZ(sizeY * sizeX),
          scale(voxelScale), invScale(1.0 / scale),
          oXc(oX + 0.5 * scale),
          oYc(oY + 0.5 * scale),
          oZc(oZ + 0.5 * scale),
          occPtr(nullptr),
          boundX(sizeX - 1),
          boundY((sizeY - 1) * stepY),
          boundZ((sizeZ - 1) * stepZ),
          stepScaleX(scale),
          stepScaleY(scale / stepY),
          stepScaleZ(scale / stepZ)
    {
        occPtr = new uint8_t[sizeVol];
        std::fill_n(occPtr, sizeVol, uint8_t(0));
    }

    ~OccupancyGridMap()
    {
        if (occPtr != nullptr)
        {
            delete[] occPtr;
        }
    }

    OccupancyGridMap(const OccupancyGridMap &) = delete;

public:
    void addWall()
    {
        for (int i = 0; i < sizeX; i++)
        {
            for (int k = 0; k < sizeZ; k++)
            {
                occPtr[i + (sizeY - 1) * stepY + k * stepZ] = 1;
                occPtr[i + 0 * stepY + k * stepZ] = 1;
            }
        }
        for (int j = 0; j < sizeY; j++)
        {
            for (int k = 0; k < sizeZ; k++)
            {
                occPtr[sizeX - 1 + j * stepY + k * stepZ] = 1;
                occPtr[0 + j * stepY + k * stepZ] = 1;
            }
        }
    }

    int getPointCloudSize()
    {
        int count = 0;
        for (int i = 0; i < sizeVol; i++)
        {
            if (occPtr[i] > 0)
                count++;
        }
        return count;
    }

    inline Eigen::Vector3i getSize(void) const
    {
        return Eigen::Vector3i(sizeX, sizeY, sizeZ);
    }

    inline double getScale(void) const
    {
        return scale;
    }

    inline Eigen::Vector3d getOrigin(void) const
    {
        return Eigen::Vector3d(oX, oY, oZ);
    }

    inline const uint8_t *getDataPtr(void) const
    {
        return occPtr;
    }

    // set the voxel to occupied, which related to pos
    inline void setOccupied(const Eigen::Vector3d &pos)
    {
        const int tempXi = (pos(0) - oX) * invScale;
        const int tempYi = (pos(1) - oY) * invScale;
        const int tempZi = (pos(2) - oZ) * invScale;
        if (tempXi >= 0 && tempYi >= 0 && tempZi >= 0 &&
            tempXi < sizeX && tempYi < sizeY && tempZi < sizeZ)
        {
            occPtr[tempXi + tempYi * stepY + tempZi * stepZ] = 1;
        }
    }

    inline void setOccupiedByEditor(const Eigen::Vector3d &pos)
    {
        const int tempXi = (pos(0) - oX) * invScale;
        const int tempYi = (pos(1) - oY) * invScale;
        const int tempZi = (pos(2) - oZ) * invScale;
        if (tempXi >= 0 && tempYi >= 0 && tempZi >= 0 &&
            tempXi < sizeX && tempYi < sizeY && tempZi < sizeZ)
        {
            occPtr[tempXi + tempYi * stepY + tempZi * stepZ] += 1;
        }
    }

    inline void setFreeByEditor(const Eigen::Vector3d &pos)
    {
        const int tempXi = (pos(0) - oX) * invScale;
        const int tempYi = (pos(1) - oY) * invScale;
        const int tempZi = (pos(2) - oZ) * invScale;
        if (tempXi >= 0 && tempYi >= 0 && tempZi >= 0 &&
            tempXi < sizeX && tempYi < sizeY && tempZi < sizeZ)
        {
            if (occPtr[tempXi + tempYi * stepY + tempZi * stepZ] > 0)
                occPtr[tempXi + tempYi * stepY + tempZi * stepZ] = uint8_t(0);
        }
    }

    inline void accumOccupied(const Eigen::Vector3d &pos)
    {
        const int tempXi = (pos(0) - oX) * invScale;
        const int tempYi = (pos(1) - oY) * invScale;
        const int tempZi = (pos(2) - oZ) * invScale;
        if (tempXi >= 0 && tempYi >= 0 && tempZi >= 0 &&
            tempXi < sizeX && tempYi < sizeY && tempZi < sizeZ)
        {
            if (occPtr[tempXi + tempYi * stepY + tempZi * stepZ] < UINT8_MAX)
            {
                occPtr[tempXi + tempYi * stepY + tempZi * stepZ]++;
            }
        }
    }

    inline void applyThresh(const int &thresh)
    {
        for (int i = 0; i < sizeVol; i++)
        {
            occPtr[i] = (occPtr[i] < thresh) ? 0 : 1;
        }
        return;
    }

    // set the voxel to occupied via 3d discrete idex
    inline void setIdxOccupied(const Eigen::Vector3i &id)
    {
        if (id(0) >= 0 && id(1) >= 0 && id(2) >= 0 &&
            id(0) < sizeX && id(1) < sizeY && id(2) < sizeZ)
        {
            occPtr[id(0) + id(1) * stepY + id(2) * stepZ] = 1;
        }
    }

    // dilate the gridmap with linear complexity
    // 线性复杂度扩张栅格地图
    inline void dilate(const double &radius)
    {
        // compute the discrete dilation radius
        int r = std::ceil(radius * invScale);
        ROS_INFO("dilate r = %d", r);
        if (r <= 0)
        {
            return;
        }
        else
        {
            // dilate the gridmap incrementally
            // store the surface in either llist or clist
            // this loop obtained the surface voxels of undilated gridmap
            std::vector<Eigen::Vector3i> llist, clist;
            llist.reserve(sizeVol);
            clist.reserve(sizeVol);

            for (int x = 0; x <= boundX; x++)
            {
                for (int y = 0; y <= boundY; y += stepY)
                {
                    for (int z = 0; z <= boundZ; z += stepZ)
                    {
                        // printf("bound Z = %d", boundZ);
                        // if (occPtr[x + y + z] == 1)
                        if (occPtr[x + y + z] == 1)
                        {
                            // this macro access 26 neibour of an occupied voxel
                            // newly occupied surface voxels are stored
                            // DILATE_LOOP(i, j, k, x, y, z,
                            //             stepY, stepZ, boundX, boundY, boundZ,
                            //             check, occPtr, idx, 2, clist)
                            my_dilate_loop(x, y, z, clist);
                        }
                    }
                }
            }
            // ROS_INFO("clist size = %d",clist.size());

            // incrementally track the obstacle surface voxels after each dilation
            for (int loop = 1; loop < r; loop++)
            {
                std::swap(clist, llist);

                for (const Eigen::Vector3i &id : llist)
                {
                    // this macro still access 26 neibour of an occupied voxel
                    // newly occupied surface voxels are stored.
                    // this newly store voxels are indeed surface of obstacles
                    // DILATE_LOOP(i, j, k, id(0), id(1), id(2),
                    //             stepY, stepZ, boundX, boundY, boundZ,
                    //             check, occPtr, idx, 2, clist)
                    my_dilate_loop(id(0), id(1), id(2), clist);
                }
            }
            dilatedSurface = clist;
        }
    }

    inline void my_dilate_loop(int x, int y, int z, std::vector<Eigen::Vector3i> &clist)
    {
        int i, j, k;
        bool check = x == 0 || x == boundX || y == 0 || y == boundY || z == 0 || z == boundZ;

        for (int ii = -1; ii <= 1; ii++)
        {
            for (int jj = -1; jj <= 1; jj++)
            {
                for (int kk = -1; kk <= 1; kk++)
                {
                    if (abs(ii) + abs(jj) + abs(kk) == 3)
                        continue;
                    i = x + ii;
                    j = y + jj * stepY;
                    k = z + kk * stepZ;
                    int offset = i + j + k;
                    if (!check && occPtr[offset] == 0)
                    {
                        // ROS_INFO("check");
                        occPtr[offset] = 2;
                        clist.emplace_back(i, j, k);
                    }
                }
            }
        }
    }

    template <typename T>
    inline void getAllOccupied(T &xyzs, const double &clipHeight) const
    {
        xyzs.clear();
        for (int x = 0; x <= boundX; x++)
        {
            for (int y = 0; y <= boundY; y += stepY)
            {
                for (int z = 0; z <= boundZ; z += stepZ)
                {
                    if (occPtr[x + y + z] && (z * stepScaleZ + oZc) < clipHeight)
                    {
                        xyzs.push_back(typename T::PointType(x * stepScaleX + oXc,
                                                             y * stepScaleY + oYc,
                                                             z * stepScaleZ + oZc));
                    }
                }
            }
        }
    }

    inline void getSurfacePointsInBox(const Eigen::Vector3i &center,
                                      const int &halfWi,
                                      std::vector<double> &xyzs) const
    {
        // this access surface points
        // it find out voxels that is closer than halfWi radius to center
        for (const Eigen::Vector3i &id : dilatedSurface)
        {
            if (abs(id(0) - center(0)) <= halfWi &&
                abs(id(1) / stepY - center(1)) <= halfWi &&
                abs(id(2) / stepZ - center(2)) <= halfWi)
            {
                xyzs.push_back(id(0) * stepScaleX + oXc);
                xyzs.push_back(id(1) * stepScaleY + oYc);
                xyzs.push_back(id(2) * stepScaleZ + oZc);
            }
        }

        return;
    }

    // query if an position is occupied by an obstacle voxel
    inline bool query(const Eigen::Vector3d &pos) const
    {
        const int tempXi = (pos(0) - oX) * invScale;
        const int tempYi = (pos(1) - oY) * invScale;
        const int tempZi = (pos(2) - oZ) * invScale;
        if (tempXi >= 0 && tempYi >= 0 && tempZi >= 0 &&
            tempXi < sizeX && tempYi < sizeY && tempZi < sizeZ)
        {
            return occPtr[tempXi + tempYi * stepY + tempZi * stepZ];
        }
        else
        {
            return true;
        }
    }

    // query if a discrete coordinate is occupied by an obstacle voxel
    inline bool queryIdx(const Eigen::Vector3i &id) const
    {
        if (id(0) >= 0 && id(1) >= 0 && id(2) >= 0 &&
            id(0) < sizeX && id(1) < sizeY && id(2) < sizeZ)
        {
            return occPtr[id(0) + id(1) * stepY + id(2) * stepZ];
        }
        else
        {
            return true;
        }
    }

    // convert the discrete coordinate to corresponding center pos
    inline Eigen::Vector3d convertPosI2D(const Eigen::Vector3i &id) const
    {
        return Eigen::Vector3d(id(0) * scale + oXc,
                               id(1) * scale + oYc,
                               id(2) * scale + oZc);
    }

    // convert the pos to its discrete coordinate
    inline Eigen::Vector3i convertPosD2I(const Eigen::Vector3d &pos) const
    {
        return Eigen::Vector3i((pos(0) - oX) * invScale,
                               (pos(1) - oY) * invScale,
                               (pos(2) - oZ) * invScale);
    }
};

// this is just a wrapper for the above class
class GlobalMap

{
public:
    Config config;
    OccupancyGridMap *ogmPtr;

    GlobalMap(const Config &conf) : config(conf), ogmPtr(nullptr) {}

    ~GlobalMap()
    {
        if (ogmPtr != nullptr)
        {
            delete ogmPtr;
        }
    }

    template <typename PointCloudType>
    inline void initialize(const PointCloudType &pc, const int &thresh)
    {
        Eigen::Vector3i xyz((config.r3Bound[1] - config.r3Bound[0]) / config.gridResolution,
                            (config.r3Bound[3] - config.r3Bound[2]) / config.gridResolution,
                            (config.r3Bound[5] - config.r3Bound[4]) / config.gridResolution);

        Eigen::Vector3d offset(config.r3Bound[0], config.r3Bound[2], config.r3Bound[4]);

        ogmPtr = new OccupancyGridMap(xyz, offset, config.gridResolution);

        size_t total = pc.size();
        // ROS_INFO("pc size = %d",total);
        double x, y, z;
        for (size_t i = 0; i < total; i++)
        {
            x = pc[i].x;
            y = pc[i].y;
            z = pc[i].z;

            if (std::isnan(x) || std::isinf(x) ||
                std::isnan(y) || std::isinf(y) ||
                std::isnan(z) || std::isinf(z))
            {
                continue;
            }
            ogmPtr->accumOccupied(Eigen::Vector3d(x, y, z));
        }
        // ROS_INFO("ogm pc size = %d",ogmPtr->getPointCloudSize());
        ogmPtr->applyThresh(thresh);
        ogmPtr->dilate(config.dilateRadius);
        ogmPtr->addWall();
        // ROS_INFO("dilated ogm pc size = %d",ogmPtr->getPointCloudSize());

        return;
    }

    template <typename PointCloudType>
    inline void initialize(const PointCloudType &pc, const int &thresh, const Eigen::Vector3d &bd_min, const Eigen::Vector3d &bd_max)
    {
        Eigen::Vector3i xyz((bd_max(0) - bd_min(0)) / config.gridResolution,
                            (bd_max(1) - bd_min(1)) / config.gridResolution,
                            (bd_max(2) - bd_min(2)) / config.gridResolution);

        Eigen::Vector3d offset(bd_min(0), bd_min(1), bd_min(2));

        ogmPtr = new OccupancyGridMap(xyz, offset, config.gridResolution);

        size_t total = pc.size();
        // ROS_INFO("pc size = %d",total);
        double x, y, z;
        for (size_t i = 0; i < total; i++)
        {
            x = pc[i].x;
            y = pc[i].y;
            z = pc[i].z;

            if (std::isnan(x) || std::isinf(x) ||
                std::isnan(y) || std::isinf(y) ||
                std::isnan(z) || std::isinf(z))
            {
                continue;
            }
            ogmPtr->accumOccupied(Eigen::Vector3d(x, y, z));
        }
        // ROS_INFO("ogm pc size = %d",ogmPtr->getPointCloudSize());
        ogmPtr->applyThresh(thresh);
        ogmPtr->dilate(config.dilateRadius);
        ogmPtr->addWall();
        // ROS_INFO("dilated ogm pc size = %d",ogmPtr->getPointCloudSize());

        return;
    }

    template <typename PointCloudType>
    inline void getPointCloud(PointCloudType &pc,
                              const double &clipHeight)
    {
        pc.clear();
        ogmPtr->getAllOccupied(pc, clipHeight);

        return;
    }

    inline bool safeQuery(const Eigen::Vector3d &p) const
    {
        return !ogmPtr->query(p);
    }

    inline Eigen::Vector3i getSize() const
    {
        return ogmPtr->getSize();
    }

    inline Eigen::Vector3d getOrigin() const
    {
        return ogmPtr->getOrigin();
    }

    inline double getScale() const
    {
        return ogmPtr->getScale();
    }
};

#endif
