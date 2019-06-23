#ifndef SFMDATA_H
#define SFMDATA_H

/**
 * @file    SFMdata.h
 * @brief   Simple example for the structure-from-motion problems
 * @author
 */

/**
 * A structure-from-motion example with landmarks
 *  - The landmarks form a 10 meter cube
 *  - The robot rotates around the landmarks, always facing towards the cube
 */

// As this is a full 3D problem, we will use Pose3 variables to represent the camera
// positions and Point3 variables (x, y, z) to represent the landmark coordinates.
// Camera observations of landmarks (i.e. pixel coordinates) will be stored as Point2 (x, y).
// We will also need a camera object to hold calibration information and perform projections.
#include "minisam/geometry/Pose3.h"
// We will also need a camera object to hold calibration information and perform projections.
#include "minisam/geometry/SimpleCamera.h"
#include "minisam/base/Matrix.h"
#include <Eigen/Core>

using namespace std;
using namespace minisam;
/* ************************************************************************* */
std::vector<Eigen::Vector3d> createPoints()
{

    // Create the set of ground-truth landmarks
    std::vector<Eigen::Vector3d> points;
    points.push_back(Eigen::Vector3d(10.0,10.0,10.0));
    points.push_back(Eigen::Vector3d(-10.0,10.0,10.0));
    points.push_back(Eigen::Vector3d(-10.0,-10.0,10.0));
    points.push_back(Eigen::Vector3d(10.0,-10.0,10.0));
    points.push_back(Eigen::Vector3d(10.0,10.0,-10.0));
    points.push_back(Eigen::Vector3d(-10.0,10.0,-10.0));
    points.push_back(Eigen::Vector3d(-10.0,-10.0,-10.0));
    points.push_back(Eigen::Vector3d(10.0,-10.0,-10.0));

    return points;
}

/* ************************************************************************* */
std::vector<Pose3> createPoses()
{

    // Create the set of ground-truth poses
    std::vector<Pose3> poses;
    double radius = 30.0;
    int i = 0;
    double theta = 0.0;
    Eigen::Vector3d up(0,0,1);
    Eigen::Vector3d target(0,0,0);
    for(; i < 8; ++i, theta += 2*M_PI/8)
    {
        Eigen::Vector3d position = Eigen::Vector3d(radius*cos(theta), radius*sin(theta), 0.0);
        SimpleCamera camera = minisam::SimpleCamera::Lookat(position, target, up);
        poses.push_back(camera.pose());
    }
    return poses;
}
/* ************************************************************************* */
#endif // SFMDATA_H
