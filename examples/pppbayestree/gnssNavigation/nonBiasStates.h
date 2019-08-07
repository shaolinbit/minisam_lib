/**
 * @file   nonBiasStates.h
 * @brief  Wrapper for all GNSS non-bias states ---> {delta pos, clock bias, clock drift, residual zenith trop}
 * @author Ryan Watson
 */

// \callgraph

#pragma once
#include <Eigen/Core>

/**
 * A vector of size 5 with multiple methods needed for GTSAM
 * @addtogroup geometry
 * \nosubgrouping
 */
class  nonBiasStates : public Eigen::VectorXd
{

private:
double x_, y_, z_, cb_, tz_;

public:

enum { dimension = 5 };

/// @name Standard Constructors
/// @{

nonBiasStates() :
        x_(0), y_(0), z_(0), cb_(0), tz_(0) {
}

/** constructor */
nonBiasStates(double x, double y, double z, double cb, double tz) :
        x_(x), y_(y), z_(z), cb_(cb), tz_(tz) {
}

// construct from 5D vector
explicit nonBiasStates(const Eigen::VectorXd& v) :
        x_(v(0)), y_(v(1)), z_(v(2)), cb_(v(3)), tz_(v(4)) {
}

// @}

/// @name Group
/// @{

/// identity for group operation
inline static nonBiasStates identity() {
        return nonBiasStates(0.0, 0.0, 0.0, 0.0, 0.0);
}

/// @}
/// @name Vector Space
/// @{

/** distance between two points
double distance(const nonBiasStates& p2, OptionalJacobian<1, 5> H1 = boost::none,
                OptionalJacobian<1, 5> H2 = boost::none) const;*/

double distance(const nonBiasStates& p2, Eigen::MatrixXd* H1 = NULL,
                Eigen::MatrixXd* H2 =NULL) const;

/** Distance of the point from the origin, with Jacobian */

double norm(Eigen::MatrixXd* H = NULL) const;


/** dot product @return this * q*/
double dot(const nonBiasStates &q, Eigen::MatrixXd* H_p = NULL,     //
           Eigen::MatrixXd* H_q = NULL) const;
/// @}
/// @name Standard Interface
/// @{

/// return as Vector5
const Eigen::VectorXd& vector() const {
        return *this;
}

/// get x
inline double x() const {
        return (*this)[0];
}

/// get y
inline double y() const {
        return (*this)[1];
}

/// get z
inline double z() const {
        return (*this)[2];
}

/// get cb
inline double cb() const {
        return (*this)[3];
}
/// @}

// get tz
inline double tz() const {
        return (*this)[4];
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
                 Eigen::MatrixXd* H1 = NULL,     //
           Eigen::MatrixXd* H2 = NULL);


/// Distance of the point from the origin, with Jacobian
double GNSS_norm5(const nonBiasStates& p, Eigen::MatrixXd* H = NULL);

/// dot product
double GNSS_dot(const nonBiasStates& p, const nonBiasStates& q,
           Eigen::MatrixXd* H_p = NULL,     //
           Eigen::MatrixXd* H_q = NULL);

template <typename A1, typename A2>
struct Range;

template <>
struct Range<nonBiasStates, nonBiasStates> {
        typedef double result_type;
        double operator()(const nonBiasStates& p, const nonBiasStates& q,
                           Eigen::MatrixXd* H1 = NULL,     //
           Eigen::MatrixXd* H2 = NULL) {
                return GNSS_distance5(p, q, H1, H2);
        }
};
