/** @file   nonBiasStates.cpp
 * @brief  non-bias state vector -- {x,y,z,cb,tz}
 * @author Ryan Watson
 */

#include "nonBiasStates.h"
#include <cmath>
#include <iostream>

using namespace std;



/* ************************************************************************* */
double nonBiasStates::distance(const nonBiasStates &q,Eigen::MatrixXd* H1,     //
           Eigen::MatrixXd* H2) const {
        return GNSS_distance5(*this,q,H1,H2);
}

double nonBiasStates::norm(Eigen::MatrixXd* H) const {
        return GNSS_norm5(*this, H);
}

double nonBiasStates::dot(const nonBiasStates &q, Eigen::MatrixXd* H1,     //
           Eigen::MatrixXd* H2) const {
        return GNSS_dot(*this, q, H1, H2);
}

/* ************************************************************************* */
ostream &operator<<(ostream &os, const nonBiasStates& p) {
        os << "   " << '[' << p.x() << ", " << p.y() << ", " << p.z() << ", " << p.cb() << ", " << p.tz() << "]'";
        return os;
}


/* ************************************************************************* */
double GNSS_distance5(const nonBiasStates &p1, const nonBiasStates &q, Eigen::MatrixXd* H1,     //
           Eigen::MatrixXd* H2) {
        double range = (q - p1).norm();
        if (H1!=NULL) {
        H1->resize(1,5);
                *H1 << p1.x() - q.x(), p1.y() - q.y(), p1.z() - q.z(), p1.cb() - q.cb(), p1.tz() - q.tz();
                *H1 = *H1 *(1. / range);
        }
        if (H2!=NULL) {
        H2->resize(1,5);
                *H2 << -p1.x() + q.x(), -p1.y() + q.y(), -p1.z() + q.z(), -q.cb() + q.cb(), -p1.tz() + q.tz();
                *H2 = *H2 *(1. / range);
        }
        return range;
}

// returns estimated pseudoragne
double GNSS_norm5(const nonBiasStates &p, Eigen::MatrixXd* H)
 {
        double r = sqrt(p.x() * p.x() + p.y() * p.y() + p.z() * p.z()) + p.cb() + p.tz();
        if (H!=NULL) {
        H->resize(1,5);
                if (fabs(r) > 1e-10)
                        *H << p.x() / r, p.y() / r, p.z() / r, 1, 1;
                else
                        *H << 1, 1, 1, 1, 1;  // really infinity, why 1 ?
        }
        return r;
}

double GNSS_dot(const nonBiasStates &p, const nonBiasStates &q, Eigen::MatrixXd* H1,     //
           Eigen::MatrixXd* H2) {
        if (H1!=NULL)
        {
        H1->resize(1,5);
        *H1 << q.x(), q.y(), q.z(), q.cb(), q.tz();
        }
        if (H2!=NULL)
        {
        H2->resize(1,5);
        *H2 << p.x(), p.y(), p.z(), p.cb(), p.tz();
        }
        return p.x() * q.x() + p.y() * q.y() + p.z() * q.z() + p.cb() * q.cb() + p.tz()*q.tz();
}

/* ************************************************************************* */
ostream &operator<<(ostream &os, const nonBiasPair &p) {
        os << p.first << " <-> " << p.second;
        return os;
}
/* ************************************************************************* */

