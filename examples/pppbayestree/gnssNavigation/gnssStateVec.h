/**
 * @file   gnssStateVec.h
 * @brief  gnss state vector ---> {delta pos, clock bias, residual zenith trop}
 * @author Ryan Watson
 */


#pragma once

#include "minisam/miniblas/minimatrix_double.h"
#include "minisam/miniblas/minivector_double.h"

/**
 * A vector of size 5 with multiple methods needed for GTSAM
 * @addtogroup geometry
 * \nosubgrouping
 */
class  gnssStateVec : public minivector
{
//public:
 //   minivector gnssstate_vector_;
//private:
 //   double x_, y_, z_, cb_, tz_;

public:

  //  enum { dimension = 5 };

/// @name Standard Constructors
/// @{

    gnssStateVec() :minivector(5)
      //  x_(0), y_(0), z_(0), cb_(0), tz_(0)
    {
       // minivector_resize(&gnssstate_vector_,5);
        minivector_set_zero(this);
    }

    /** constructor */
    gnssStateVec(double x, double y, double z, double cb, double tz) :minivector(5)
       // x_(x), y_(y), z_(z), cb_(cb), tz_(tz)
    {
        //gnssstate_vector_=NULL_VECTOR;
        //minivector_resize(&gnssstate_vector_,5);
        data[0]=x;
        data[1]=y;
        data[2]=z;
        data[3]=cb;
        data[4]=tz;
    }

/// construct from 3D vector
    explicit gnssStateVec(const minivector& v) :minivector(v)
        //x_(v.data[0]), y_(v.data[1]), z_(v.data[2]), cb_(v.data[3]), tz_(v.data[4])
    {
       // minivector_resize(&gnssstate_vector_,5);
       // minivector_memcpy(&gnssstate_vector_,v);
    }
    gnssStateVec(const minivector* v) :minivector(v)
    {

    }


/// @}

/// @name Group
/// @{

/// identity for group operation
    inline static gnssStateVec identity()
    {
        return gnssStateVec(0.0, 0.0, 0.0, 0.0, 0.0);
    }

/// @}
/// @name Vector Space
/// @{


    /** distance between two points */
    double distance(const gnssStateVec& p2, minimatrix* H1 = NULL,
                    minimatrix*  H2 = NULL) const;

    /** Distance of the point from the origin, with Jacobian */
//double norm(OptionalJacobian<1,5> H = boost::none) const;
    /** Distance of the point from the origin, with Jacobian */
    double norm(minimatrix* H1 = NULL) const;


    double dot(const gnssStateVec& q, minimatrix* H_p= NULL,
               minimatrix*  H_q = NULL) const;
/// @}
/// @name Standard Interface
/// @{

/// return as Vector5
    const minivector& vector() const
    {
        return *this;
    }

/// get x
    inline double x() const
    {
        return minivector_get(this,0);
    }

/// get y
    inline double y() const
    {
        return  minivector_get(this,1);
    }

/// get z
    inline double z() const
    {
        return  minivector_get(this,2);
    }

/// get cb
    inline double cb() const
    {
        return  minivector_get(this,3);
    }
/// @}

// get tz
    inline double tz() const
    {
        return  minivector_get(this,4);
    }
/// @}

    void setgnssStateVec(double x, double y, double z, double cb, double tz)
    {
        /*
        x_=x;
        y_=y;
        z_=z;
        cb_=cb;
        tz_=tz;
        minivector_resize(&gnssstate_vector_,5);

        gnssstate_vector_.data[0]=x_;
        gnssstate_vector_.data[1]=y_;
        gnssstate_vector_.data[2]=z_;
        gnssstate_vector_.data[3]=cb_;
        gnssstate_vector_.data[4]=tz_;*/
        minivector_set(this,0,x); minivector_set(this,1,y); minivector_set(this,2,z);
        minivector_set(this,3,cb); minivector_set(this,4,tz);

    }

    void setgnssStateVec(const minivector& v)
    {
        /*
        x_=v.data[0];
        y_=v.data[1];
        z_=v.data[2];
        cb_=v.data[3];
        tz_=v.data[4];
        minivector_resize(&gnssstate_vector_,5);
        gnssstate_vector_.data[0]=x_;
        gnssstate_vector_.data[1]=y_;
        gnssstate_vector_.data[2]=z_;
        gnssstate_vector_.data[3]=cb_;
        gnssstate_vector_.data[4]=tz_;
        */
        minivector_memcpy(this,v);
    }

/// Output stream operator
    friend std::ostream &operator<<(std::ostream &os, const gnssStateVec& p);



};


// Convenience typedef
typedef std::pair<gnssStateVec, gnssStateVec> Point5Pair;
std::ostream &operator<<(std::ostream &os, const Point5Pair &p);

/// distance between two points

double GNSS_distance5(const gnssStateVec& p1, const gnssStateVec& q,
                      minimatrix* H1 = NULL,
                      minimatrix* H2 = NULL);


/// Distance of the point from the origin, with Jacobian

double GNSS_norm5(const gnssStateVec& p, minimatrix* H = NULL);

/// dot product

double GNSS_dot(const gnssStateVec& p, const gnssStateVec& q,
                minimatrix* H_p= NULL,
                minimatrix*  H_q = NULL);

template <typename A1, typename A2>
struct Range;

template <>
struct Range<gnssStateVec, gnssStateVec>
{
    typedef double result_type;

    double operator()(const gnssStateVec& p, const gnssStateVec& q,
                      minimatrix* H1 = NULL,
                      minimatrix* H2 = NULL)
    {
        return GNSS_distance5(p, q, H1, H2);
    }
};
