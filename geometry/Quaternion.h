/**
 * @file   Quaternion.h
 * @brief  Lie Group wrapper for Eigen Quaternions
 * @author
 **/

#pragma once

#include "../geometry/SO3.h" // Logmap/Expmap derivatives
#include "../base/Matrix.h"
#include <limits>
#include <iostream>

class QQuaternion
{
    /// @name Group traits
    /// @{
public:
    Eigen::Quaternion<double, Eigen::DontAlign> Q_;
public:
    QQuaternion(double q0,double q1,double q2,double q3)
    {
        Q_=Eigen::Quaternion<double, Eigen::DontAlign>(q0,q1,q2,q3);
    }
    QQuaternion(const Eigen::Matrix3d& df)
    {
        Q_=Eigen::Quaternion<double, Eigen::DontAlign>(df);
    }
     QQuaternion(const Eigen::Quaternion<double, Eigen::DontAlign>& df)
    {
        Q_=df;
    }
      QQuaternion(const Eigen::AngleAxisd& df)
    {
        Q_=df;
    }


    QQuaternion() {}
    ~QQuaternion() {}
public:
    static QQuaternion Identity()
    {
        QQuaternion qq;
        qq.Q_=Eigen::Quaternion<double, Eigen::DontAlign>::Identity();
        return qq;
    }

    /// @}
    /// @name Basic manifold traits
    /// @{
    enum
    {
        dimension = 3
    };
    //typedef OptionalJacobian<3, 3> ChartJacobian;
    //typedef Eigen::Matrix<_Scalar, 3, 1, _Options, 3, 1> TangentVector;

    /// @}
    /// @name Lie group traits
    /// @{
    static QQuaternion Compose(const QQuaternion &g, const QQuaternion & h)
    {
        // if (Hg) *Hg = h.toRotationMatrix().transpose();
        // if (Hh) *Hh = I_3x3;
        //return g * h;
        QQuaternion qq;
        qq.Q_=g.Q_*h.Q_;
        return qq;
    }
    static QQuaternion Compose(const QQuaternion &g, const QQuaternion & h,
                               Eigen::Matrix3d* Hg, Eigen::Matrix3d* Hh)
    {
        //if (Hg)
        *Hg = h.Q_.toRotationMatrix().transpose();
        // if (Hh)
        *Hh = Eigen::Matrix3d::Identity();
        QQuaternion qq;
        qq.Q_=g.Q_*h.Q_;
        return qq;
        //  return g * h;
    }

    static QQuaternion Between(const QQuaternion &g, const QQuaternion & h)
    {
        QQuaternion d;
        d.Q_ = g.Q_.inverse() * h.Q_;
        // if (Hg) *Hg = -d.toRotationMatrix().transpose();
        // if (Hh) *Hh = I_3x3;
        return d;
    }

    static QQuaternion Between(const QQuaternion &g, const QQuaternion & h,
                               Eigen::Matrix3d* Hg, Eigen::Matrix3d* Hh)
    {
        QQuaternion d;
        d.Q_ = g.Q_.inverse() * h.Q_;
        //if (Hg)
        *Hg = -d.Q_.toRotationMatrix().transpose();
        // if (Hh)
        *Hh = Eigen::Matrix3d::Identity();
        return d;
    }

    static QQuaternion Inverse(const QQuaternion &g)
    {
        // if (H) *H = -g.toRotationMatrix();
        QQuaternion dg;
        dg.Q_=g.Q_.inverse();
        return dg;
        // return g.inverse();
    }

    static QQuaternion Inverse(const QQuaternion &g,
                               Eigen::Matrix3d* H)
    {
        //if (H)
        *H = -g.Q_.toRotationMatrix();
        QQuaternion dg;
        dg.Q_=g.Q_.inverse();
        return dg;
    }

    /// Exponential map, using the inlined code from Eigen's conversion from axis/angle
    static QQuaternion Expmap(const Eigen::Ref<const Eigen::Vector3d>& omega)
    {
        using std::cos;
        using std::sin;
        //if (H) *H = SO3::ExpmapDerivative(omega.template cast<double>());
        double theta2 = omega.dot(omega);
        if (theta2 > std::numeric_limits<double>::epsilon())
        {
            double theta = std::sqrt(theta2);
            double ha = 0.5 * theta;
            Eigen::Vector3d vec = (sin(ha) / theta) * omega;
            return QQuaternion(cos(ha), vec.x(), vec.y(), vec.z());
        }
        else
        {
            // first order approximation sin(theta/2)/theta = 0.5
            Eigen::Vector3d vec = 0.5 * omega;
            return QQuaternion(1.0, vec.x(), vec.y(), vec.z());
        }
    }

    static QQuaternion Expmap(const Eigen::Ref<const Eigen::Vector3d>& omega,
                              Eigen::Matrix3d* H)
    {
        using std::cos;
        using std::sin;
        *H = SO3::ExpmapDerivative(omega.template cast<double>());
        double theta2 = omega.dot(omega);
        if (theta2 > std::numeric_limits<double>::epsilon())
        {
            double theta = std::sqrt(theta2);
            double ha = 0.5 * theta;
            Eigen::Vector3d vec = (sin(ha) / theta) * omega;
            return QQuaternion(cos(ha), vec.x(), vec.y(), vec.z());
        }
        else
        {
            // first order approximation sin(theta/2)/theta = 0.5
            Eigen::Vector3d vec = 0.5* omega;
            return QQuaternion(1.0, vec.x(), vec.y(), vec.z());
        }
    }

    /// We use our own Logmap, as there is a slight bug in Eigen
    static Eigen::Vector3d Logmap(const QQuaternion& q)
    {
        using std::acos;
        using std::sqrt;

        // define these compile time constants to avoid std::abs:
        static const double twoPi = 2.0 * M_PI, NearlyOne = 1.0 - 1e-10,
                            NearlyNegativeOne = -1.0 + 1e-10;

        Eigen::Vector3d omega;

        const double qw = q.Q_.w();
        // See Quaternion-Logmap.nb in doc for Taylor expansions
        if (qw > NearlyOne)
        {
            // Taylor expansion of (angle / s) at 1
            // (2 + 2 * (1-qw) / 3) * q.vec();
            omega = ( 8. / 3. - 2. / 3. * qw) * q.Q_.vec();
        }
        else if (qw < NearlyNegativeOne)
        {
            // Taylor expansion of (angle / s) at -1
            // (-2 - 2 * (1 + qw) / 3) * q.vec();
            omega = (-8. / 3. - 2. / 3. * qw) * q.Q_.vec();
        }
        else
        {
            // Normal, away from zero case
            double angle = 2 * acos(qw), s = sqrt(1 - qw * qw);
            // Important:  convert to [-pi,pi] to keep error continuous
            if (angle > M_PI)
                angle -= twoPi;
            else if (angle < -M_PI)
                angle += twoPi;
            omega = (angle / s) * q.Q_.vec();
        }

        // if(H) *H = SO3::LogmapDerivative(omega.template cast<double>());
        return omega;
    }

    static Eigen::Vector3d Logmap(const QQuaternion& q, Eigen::Matrix3d* H)
    {
        using std::acos;
        using std::sqrt;

        // define these compile time constants to avoid std::abs:
        static const double twoPi = 2.0 * M_PI, NearlyOne = 1.0 - 1e-10,
                            NearlyNegativeOne = -1.0 + 1e-10;

        Eigen::Vector3d omega;

        const double qw = q.Q_.w();
        // See Quaternion-Logmap.nb in doc for Taylor expansions
        if (qw > NearlyOne)
        {
            // Taylor expansion of (angle / s) at 1
            // (2 + 2 * (1-qw) / 3) * q.vec();
            omega = ( 8. / 3. - 2. / 3. * qw) * q.Q_.vec();
        }
        else if (qw < NearlyNegativeOne)
        {
            // Taylor expansion of (angle / s) at -1
            // (-2 - 2 * (1 + qw) / 3) * q.vec();
            omega = (-8. / 3. - 2. / 3. * qw) * q.Q_.vec();
        }
        else
        {
            // Normal, away from zero case
            double angle = 2 * acos(qw), s = sqrt(1 - qw * qw);
            // Important:  convert to [-pi,pi] to keep error continuous
            if (angle > M_PI)
                angle -= twoPi;
            else if (angle < -M_PI)
                angle += twoPi;
            omega = (angle / s) * q.Q_.vec();
        }

        //if(H)
        *H = SO3::LogmapDerivative(omega.template cast<double>());
        return omega;
    }

    /// @}
    /// @name Manifold traits
    /// @{

    static Eigen::Vector3d Local(const QQuaternion& g, const QQuaternion& h)
    {
        QQuaternion b = Between(g, h);
        Eigen::Matrix3d D_v_b;
        Eigen::Vector3d v = Logmap(b);
        //if (H1) *H1 = D_v_b * (*H1);
        // if (H2) *H2 = D_v_b * (*H2);
        return v;
    }

    static Eigen::Vector3d Local(const QQuaternion& g, const QQuaternion& h,
                                 Eigen::Matrix3d* H1, Eigen::Matrix3d* H2)
    {
        QQuaternion b = Between(g, h, H1, H2);
        Eigen::Matrix3d D_v_b;
        Eigen::Vector3d v = Logmap(b, &D_v_b);
        // if (H1)
        *H1 = D_v_b * (*H1);
        //if (H2)
        *H2 = D_v_b * (*H2);
        return v;
    }

    static QQuaternion Retract(const QQuaternion& g, const Eigen::Vector3d& v)
    {
        Eigen::Matrix3d D_h_v;
        QQuaternion b = Expmap(v);
        QQuaternion h = Compose(g, b);
        // if (H2) *H2 = (*H2) * D_h_v;
        return h;
    }
    static QQuaternion Retract(const QQuaternion& g, const Eigen::Vector3d& v,
                               Eigen::Matrix3d* H1, Eigen::Matrix3d* H2)
    {
        Eigen::Matrix3d D_h_v;
        QQuaternion b = Expmap(v,&D_h_v);
        QQuaternion h = Compose(g, b, H1, H2);
        //if (H2)
        *H2 = (*H2) * D_h_v;
        return h;
    }

    inline QQuaternion& operator=(const QQuaternion& rObj)
    {
        Q_=rObj.Q_;
    }

};


