/*
 * @file OrientedPlane3.h
 * @brief An infinite plane, represented by a normal direction and perpendicular distance
 */

#pragma once

#include "../geometry/Rot3.h"
#include "../geometry/Unit3.h"
#include "../geometry/Pose3.h"

namespace minisam
{

/**
 * @brief Represents an infinite plane in 3D, which is composed of a planar normal and its
 *  perpendicular distance to the origin.
 * Currently it provides a transform of the plane, and a norm 1 differencing of two planes.
 * Refer to Trevor12iros for more math details.
 */
class  OrientedPlane3:public minivector
{
public:

    /// @name Constructors
    /// @{

    /// Default constructor
    OrientedPlane3() :minivector(4)
    {
        data[0]=1.0;
        data[1]=0.0;
        data[2]=0.0;
        data[3]=0.0;
        dimension = 3;
    }

    /// Construct from a Unit3 and a distance
    OrientedPlane3(const Unit3& s, double d) :minivector(4)
    {
        data[0]=s.data[0];
        data[1]=s.data[1];
        data[2]=s.data[2];
        data[3]=s.data[3];
        dimension = 3;
    }

    /// Construct from a vector of plane coefficients
    OrientedPlane3(const minivector& vec) :minivector(vec)
    {
        dimension = 3;
    }

    /// Construct from four numbers of plane coeffcients (a, b, c, d)
    OrientedPlane3(double a, double b, double c, double d):minivector(4)
    {
        data[0]=a;
        data[1]=b;
        data[2]=c;
        data[3]=d;
        dimension = 3;
    }

    /// @}

    /** Transforms a plane to the specified pose
     * @param xr a transformation in current coordiante
     * @param Hp optional Jacobian wrpt the destination plane
     * @param Hr optional jacobian wrpt the pose transformation
     * @return the transformed plane
     */
    OrientedPlane3 transform(const Pose3& xr,minimatrix* Hp=NULL,//3*3
                             minimatrix* Hr=NULL//3*6
                            ) const;

    /**
     * @deprecated the static method has wrong Jacobian order,
     *    please use the member method transform()
     * @param The raw plane
     * @param xr a transformation in current coordiante
     * @param Hr optional jacobian wrpt the pose transformation
     * @param Hp optional Jacobian wrpt the destination plane
     * @return the transformed plane
     */
    static OrientedPlane3 Transform(const OrientedPlane3& plane,
                                    const Pose3& xr,minimatrix* Hp=NULL,//3*3
                                    minimatrix* Hr=NULL//3*6
                                   )
    {
        return plane.transform(xr, Hp, Hr);
    }

    /** Computes the error between two planes.
     *  The error is a norm 1 difference in tangent space.
     * @param the other plane
     */
    minivector error(const OrientedPlane3& plane) const;


    /** Computes the error between the two planes, with derivatives.
     *  This uses Unit3::errorVector, as opposed to the other .error() in this class, which uses
     *  Unit3::localCoordinates. This one has correct derivatives.
     *  NOTE(hayk): The derivatives are zero when normals are exactly orthogonal.
     * @param the other plane
     */
    minivector errorVector(const OrientedPlane3& other, minimatrix* H1 = NULL, //3*3
                           minimatrix* H2 = NULL //3*3
                          ) const;

    /// Dimensionality of tangent space = 3 DOF
    inline size_t dim() const
    {
        return 3;
    }

    /// The retract function
    virtual minimatrix* Retract(const minimatrix* mpose);

    /// The local coordinates function
    virtual minimatrix LocalCoordinates(const minimatrix* mpose) const;

    /// Returns the plane coefficients
    inline minivector planeCoefficients() const
    {
        return *this;
    }

    /// Return the normal
    inline Unit3 normal() const
    {
        return Unit3(data[0],data[1],data[2]);
    }

    /// Return the perpendicular distance to the origin
    inline double distance() const
    {
        return data[3];
    }

    /// @}
};


} // namespace gtsam

