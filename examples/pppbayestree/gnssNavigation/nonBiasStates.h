/**
 * @file   nonBiasStates.h
 * @brief  Wrapper for all GNSS non-bias states ---> {delta pos, clock bias, clock drift, residual zenith trop}
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
class  nonBiasStates:public minivector
{


public:


/// @name Standard Constructors
/// @{

    nonBiasStates() :minivector(5,0.0)
    {
       
    }

    /** constructor */
    nonBiasStates(double x, double y, double z, double cb, double tz) :minivector(5)
    {
       
        data[0]=x;
        data[1]=y;
        data[2]=z;
        data[3]=cb;
        data[4]=tz;

    }

// construct from 5D vector
    explicit nonBiasStates(const minivector& v) :minivector(v)
    {

    }

    ~nonBiasStates()
    {
       
    }

// @}

/// @name Group
/// @{

/// identity for group operation
    inline static nonBiasStates identity()
    {
        return nonBiasStates(0.0, 0.0, 0.0, 0.0, 0.0);
    }

/// @}
/// @name Vector Space
/// @{

    double distance(const nonBiasStates& p2, minimatrix* H1 = NULL,
                    minimatrix* H2 =NULL) const;

    /** Distance of the point from the origin, with Jacobian */

    double norm(minimatrix* H = NULL) const;


    /** dot product @return this * q*/
    double dot(const nonBiasStates &q, minimatrix* H_p = NULL,     //
               minimatrix* H_q = NULL) const;
/// @}
/// @name Standard Interface
/// @{

/// return as Vector5
    minivector vector() const
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
        return minivector_get(this,1);
    }

/// get z
    inline double z() const
    {
        return minivector_get(this,2);
    }

/// get cb
    inline double cb() const
    {
        return minivector_get(this,3);
    }
/// @}

// get tz
    inline double tz() const
    {
        return minivector_get(this,4);
    }
/// @}

/// Output stream operator
    friend std::ostream &operator<<(std::ostream &os, const nonBiasStates& p);

};

// Convenience typedef
typedef std::pair<nonBiasStates, nonBiasStates> nonBiasPair;
std::ostream &operator<<(std::ostream &os, const nonBiasPair &p);

/// distance between two points
double GNSS_distance5(const nonBiasStates& p1, const nonBiasStates& q,
                      minimatrix* H1 = NULL,     //
                      minimatrix* H2 = NULL);


/// Distance of the point from the origin, with Jacobian
double GNSS_norm5(const nonBiasStates& p, minimatrix* H = NULL);

/// dot product
double GNSS_dot(const nonBiasStates& p, const nonBiasStates& q,
                minimatrix* H_p = NULL,     //
                minimatrix* H_q = NULL);

template <typename A1, typename A2>
struct Range;

template <>
struct Range<nonBiasStates, nonBiasStates>
{
    typedef double result_type;
    double operator()(const nonBiasStates& p, const nonBiasStates& q,
                      minimatrix* H1 = NULL,     //
                      minimatrix* H2 = NULL)
    {
        return GNSS_distance5(p, q, H1, H2);
    }
};
