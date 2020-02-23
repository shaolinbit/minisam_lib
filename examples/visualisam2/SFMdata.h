#ifndef SFMDATA_H
#define SFMDATA_H

/**
 * @file    SFMdata.h
 * @brief   Simple example for the structure-from-motion problems
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
#include "geometry/Pose3.h"

// We will also need a camera object to hold calibration information and perform projections.
#include "geometry/SimpleCamera.h"
#include "mat/Matrix.h"

using namespace std;
using namespace minisam;
/* ************************************************************************* */
std::vector<minivector*> createPoints()
{

    // Create the set of ground-truth landmarks
    std::vector<minivector*> points;
    points.push_back(new minivector(10.0,10.0,10.0));
    points.push_back(new minivector(-10.0,10.0,10.0));
    points.push_back(new minivector(-10.0,-10.0,10.0));
    points.push_back(new minivector(10.0,-10.0,10.0));
    points.push_back(new minivector(10.0,10.0,-10.0));
    points.push_back(new minivector(-10.0,10.0,-10.0));
    points.push_back(new minivector(-10.0,-10.0,-10.0));
    points.push_back(new minivector(10.0,-10.0,-10.0));

    return points;
}
void SFMData_ClearPoints(std::vector<minivector*>& cpoints)
{
   for(minivector* cpppoint:cpoints)
   {
     if(cpppoint!=NULL)
        {
        delete cpppoint;
        cpppoint=NULL;
        }
   }
   cpoints.clear();

}

/* ************************************************************************* */
std::vector<Pose3> createPoses()
{

    // Create the set of ground-truth poses
    std::vector<Pose3> poses;
    double radius = 30.0;
    int i = 0;
    double theta = 0.0;
    minivector up(0,0,1);
    minivector target(0,0,0);
    for(; i < 8; ++i, theta += 2*M_PI/8)
    {
        minivector position(radius*cos(theta), radius*sin(theta), 0.0);
       // SimpleCamera camera = SimpleCamera::Lookat(position, target, up);
        SimpleCamera camera(PinholeBaseLookatPose(position, target, up), Cal3_S2());
        poses.push_back(camera.pose());
    }
    return poses;
}
/* ************************************************************************* */
#endif // SFMDATA_H
