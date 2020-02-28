
#pragma once

#include "../mat/MatCal.h"
#include "../mat/MatCal.h"
#include <string>

namespace minisam
{

/// Represents a 3D point on a unit sphere.
class Unit3:public minivector
{

public:


    /// @name Constructors
    /// @{

    /// Default constructor
    Unit3():minivector(3)
    {
        data[0]=1.0;
        data[1]=0.0;
        data[2]=0.0;
        dimension=2;
    }
    ~Unit3()
    {

    }

    /// Construct from point
    explicit Unit3(const minivector& p):minivector(3)
    {
        normalize3d(p,this);
        dimension=2;

    }

    /// Construct from x,y,z
    Unit3(double x, double y, double z):minivector(3)
    {
        data[0]=x;
        data[1]=y;
        data[2]=z;
        dimension=2;
        normalize3d(this);

    }

    /// Construct from 2D point in plane at focal length f
    /// Unit3(p,1) can be viewed as normalized homogeneous coordinates of 2D point
    explicit Unit3(const minivector& p, double f):minivector(3)
    {
        data[0]=p.data[0];
        data[1]=p.data[p.prd];
        data[2]=f;
        normalize3d(this);
    }

    /// Copy constructor
    Unit3(const Unit3& u):minivector(3)
    {
        minivector_memcpy(this,u);
    }

    /// Copy assignment
    Unit3& operator=(const Unit3 & u)
    {
        minivector_memcpy(this,u);
        return *this;
    }

    /// Named constructor from Point3 with optional Jacobian
    static Unit3 FromPoint3(const minivector& point);
    static Unit3 FromPoint3(const minivector& point, //
                            minimatrix* H);

    /// @}

    /// @name Testable
    /// @{

    friend std::ostream& operator<<(std::ostream& os, const Unit3& pair);

    /// @}

    /// @name Other functionality
    /// @{

    /**
     * Returns the local coordinate frame to tangent plane
     * It is a 3*2 matrix [b1 b2] composed of two orthogonal directions
     * tangent to the sphere at the current direction.
     * Provides derivatives of the basis with the two basis vectors stacked up as a 6x1.
     */
     minimatrix basis() const;
     minimatrix basis(minimatrix* H) ;
    /// Return skew-symmetric associated with 3D point on unit sphere
    minimatrix skew() const;

    /// Return unit-norm Point3
    minivector point3() const;
    minivector point3(minimatrix* H) ;
    /// Return unit-norm Vector
    minivector unitVector() const;
    minivector unitVector(minimatrix* H) ;
    /// Return scaled direction as Point3
    friend minivector operator*(double s, const Unit3& d)
    {
        minivector p3(3);
        minivector_scale_vec(&p3,s,d);

        return p3;
    }

    /// Return dot product with q
    double dot(const Unit3& q) const;
    double dot( Unit3& q, minimatrix* H1, //
                minimatrix* H2) ;

    /// Signed, vector-valued error between two directions
    /// @deprecated, errorVector has the proper derivatives, this confusingly has only the second.
    minivector error(const Unit3& q) ;
    minivector error(Unit3& q, minimatrix* H_q);
    /// Signed, vector-valued error between two directions
    /// NOTE(hayk): This method has zero derivatives if this (p) and q are orthogonal.
    minivector errorVector(const Unit3& q);

    minivector errorVector( Unit3& q, minimatrix* H_p, //
                            minimatrix* H_q) ;

    /// Distance between two directions
    double distance(const Unit3& q) ;
    double distance( Unit3& q, minimatrix* H) ;
    /// Cross-product between two Unit3s
    Unit3 cross(const Unit3& q) const
    {
        minivector bp=cross3d(*this,q);
        Unit3 result(bp);
        return result;
    }

    /// Cross-product w Point3
    minivector cross(const minivector& q) const
    {
        return cross3d(point3(),q);
    }

    /// @}

    /// @name Manifold
    /// @{


    /// Dimensionality of tangent space = 2 DOF
    inline int dim() const
    {
        return 2;
    }

    enum CoordinatesMode
    {
        EXPMAP, ///< Use the exponential map to retract
        RENORM ///< Retract with vector addition and renormalize.
    };


    /// The retract function
    virtual minimatrix* Retract(const minimatrix* mpose);

    /// The local coordinates function
    virtual minimatrix LocalCoordinates(const minimatrix* mpose,minimatrix* H1=NULL,minimatrix* H2=NULL) const;

    /// @}

};
};
