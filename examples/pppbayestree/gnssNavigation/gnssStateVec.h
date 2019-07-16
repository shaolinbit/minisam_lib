/**
 * @file   gnssStateVec.h
 * @brief  gnss state vector ---> {delta pos, clock bias, residual zenith trop}
 * @author Ryan Watson
 */

// \callgraph

#pragma once

#include <Eigen/Core>


/**
 * A vector of size 5 with multiple methods
 * @addtogroup geometry
 * \nosubgrouping
 */
class  gnssStateVec : public Eigen::VectorXd
 {

private:
double x_, y_, z_, cb_, tz_;

public:

enum { dimension = 5 };

/// @name Standard Constructors
/// @{

gnssStateVec() :
        x_(0), y_(0), z_(0), cb_(0), tz_(0) {
}

/** constructor */
gnssStateVec(double x, double y, double z, double cb, double tz) :
        x_(x), y_(y), z_(z), cb_(cb), tz_(tz) {
}

/// construct from 3D vector
explicit gnssStateVec(const Eigen::VectorXd& v) :
        x_(v(0)), y_(v(1)), z_(v(2)), cb_(v(3)), tz_(v(4)) {
}


/// @}

/// @name Group
/// @{

/// identity for group operation
inline static gnssStateVec identity() {
        return gnssStateVec(0.0, 0.0, 0.0, 0.0, 0.0);
}

/// @}
/// @name Vector Space
/// @{


/** distance between two points */
double distance(const gnssStateVec& p2, Eigen::MatrixXd* H1 = NULL,
                Eigen::MatrixXd*  H2 = NULL) const;

/** Distance of the point from the origin, with Jacobian */
//double norm(OptionalJacobian<1,5> H = boost::none) const;
/** Distance of the point from the origin, with Jacobian */
double norm(Eigen::MatrixXd* H1 = NULL) const;


double dot(const gnssStateVec& q, Eigen::MatrixXd* H_p= NULL,
                Eigen::MatrixXd*  H_q = NULL) const;
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
friend std::ostream &operator<<(std::ostream &os, const gnssStateVec& p);



};


// Convenience typedef
typedef std::pair<gnssStateVec, gnssStateVec> Point5Pair;
std::ostream &operator<<(std::ostream &os, const Point5Pair &p);

/// distance between two points

double GNSS_distance5(const gnssStateVec& p1, const gnssStateVec& q,
                Eigen::MatrixXd* H1 = NULL,
                Eigen::MatrixXd* H2 = NULL);


/// Distance of the point from the origin, with Jacobian

double GNSS_norm5(const gnssStateVec& p, Eigen::MatrixXd* H = NULL);

/// dot product

double GNSS_dot(const gnssStateVec& p, const gnssStateVec& q,
           Eigen::MatrixXd* H_p= NULL,
          Eigen::MatrixXd*  H_q = NULL);

template <typename A1, typename A2>
struct Range;

template <>
struct Range<gnssStateVec, gnssStateVec> {
        typedef double result_type;

    double operator()(const gnssStateVec& p, const gnssStateVec& q,
                           Eigen::MatrixXd* H1 = NULL,
                Eigen::MatrixXd* H2 = NULL) {
                return GNSS_distance5(p, q, H1, H2);
        }
};

