#ifndef QUATERNION_H
#define QUATERNION_H


/**
 * @file   Quaternion4.h
 * @brief  Lie Group wrapper for Quaternion4s
 **/


#pragma once

#include "../geometry/SO3.h"
#include "../mat/Matrix.h"
#include "../mat/MatCal.h"
#include "../miniblas/minivector_double.h"
#include <limits>
#include <iostream>
namespace minisam
{

class Quaternion4:public minivector
{
public:
    Quaternion4(double q0,double q1,double q2,double q3):minivector(4)
    {
        dimension=3;
        data[0]=q0;
        data[1]=q1;
        data[2]=q2;
        data[3]=q3;
    }
    Quaternion4(const minimatrix& DCM):minivector(4)
    {
        double q0, q1, q2, q3, qq4;
        dimension=3;
        if(DCM.data[0]>=DCM.data[4]+DCM.data[8])
        {
            q1 =0.5*sqrt(1+DCM.data[0]-DCM.data[4]-DCM.data[8]);
            qq4 = 0.25/q1;
            q0 = (DCM.data[6]-DCM.data[5])*qq4;
            q2 = (DCM.data[1]+DCM.data[3])*qq4;
            q3 =(DCM.data[2]+DCM.data[6])*qq4;
        }
        else if(DCM.data[4]>=DCM.data[0]+DCM.data[8])
        {
            q2 =0.5*sqrt(1-DCM.data[0]+DCM.data[4]-DCM.data[8]);
            qq4 = 0.25/q2;
            q0 = (DCM.data[2]-DCM.data[6])*qq4;
            q1 = (DCM.data[1]+DCM.data[3])*qq4;
            q3 = (DCM.data[5]+DCM.data[7])*qq4;
        }
        else if(DCM.data[8]>=DCM.data[0]+DCM.data[4])
        {
            q3 =0.5*sqrt(1-DCM.data[0]-DCM.data[4]+DCM.data[8]);
            qq4 = 0.25/q3;
            q0 = (DCM.data[3]-DCM.data[1])*qq4;
            q1 = (DCM.data[2]+DCM.data[6])*qq4;
            q2 = (DCM.data[5]+DCM.data[7])*qq4;
        }
        else
        {
            q0 =0.5*sqrt(1+DCM.data[0]+DCM.data[4]+DCM.data[8]);
            qq4 = 0.25/q0;
            q1 = (DCM.data[7]-DCM.data[5])*qq4;
            q2 = (DCM.data[2]-DCM.data[6])*qq4;
            q3 = (DCM.data[3]-DCM.data[1])*qq4;
        }
        double nq = 1/sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
        data[0]=q0 * nq;
        data[1]=q1 * nq;
        data[2]=q2 *nq;
        data[3]= q3 * nq;
    }
    Quaternion4(const minivector& df):minivector(df)
    {
        dimension=3;
    }

    Quaternion4(const Quaternion4& qq):minivector(4)
    {
        dimension=3;
        minivector_memcpy(this,qq);
    }

    Quaternion4():minivector(4)
    {

        data[0]=1.0;
        data[1]=0.0;
        data[2]=0.0;
        data[3]=0.0;
        dimension=3;
    }
    ~Quaternion4()
    {

    }
public:
    static Quaternion4 Identity()
    {
        Quaternion4 qq;
        return qq;
    }

    /// @}
    /// @name Basic manifold traits
    /// @{


    /// @}
    /// @name Lie group traits
    /// @{
    static Quaternion4 Compose(const Quaternion4 &g, const Quaternion4 & h)
    {
        Quaternion4 qq;
        Quaternion_Multiply(&qq,g,h);
        return qq;
    }
    static Quaternion4 Compose(const Quaternion4 &g, const Quaternion4 & h,
                               minimatrix* Hg, minimatrix* Hh)
    {
        QuaternionToMatrixTrans(h,Hg);
        minimatrix_resize(Hh,3,3);
        minimatrix_set_identity(Hh);

        Quaternion4 qq;
        Quaternion_Multiply(&qq,g,h);
        return qq;
    }

    virtual minimatrix between(const minimatrix* mpose) const
    {
        Quaternion4 d;
        Quaternion_inv_Multiply(&d,*this,*mpose);
        return d;

    }
    virtual minimatrix between(const minimatrix* mpose,minimatrix& H1,minimatrix& H2) const
    {
        Quaternion4 d;
        Quaternion_inv_Multiply(&d,*this,*mpose);

        minimatrix_resize(&H1,3,3);
        QuaternionToMatrixTrans(d,&H1);
        minimatrix_scale(&H1,-1);


        minimatrix_resize(&H2,3,3);
        minimatrix_set_identity(&H2);


        return d;

    }

    static Quaternion4 Between(const Quaternion4 &g, const Quaternion4 & h)
    {
        Quaternion4 d;
        Quaternion_inv_Multiply(&d,g,h);
        return d;
    }

    static Quaternion4 Between(const Quaternion4 &g, const Quaternion4 & h,
                               minimatrix* Hg, minimatrix* Hh)
    {
        Quaternion4 d;
        Quaternion_inv_Multiply(&d,g,h);

        minimatrix_resize(Hg,3,3);
        QuaternionToMatrixTrans(d,Hg);
        minimatrix_scale(Hg,-1);

        minimatrix_resize(Hh,3,3);
        minimatrix_set_identity(Hh);
        return d;
    }

    static Quaternion4 Inverse(const Quaternion4 &g)
    {
        Quaternion4 dg;
        dg.data[0]=g.data[g.prd*0];
        dg.data[1]=-g.data[g.prd*1];
        dg.data[2]=-g.data[g.prd*2];
        dg.data[3]=-g.data[g.prd*3];

        return dg;
    }

    static Quaternion4 Inverse(const Quaternion4 &g,
                               minimatrix* H)
    {
        QuaternionToMatrix(g,H);
        minimatrix_scale(H,-1.0);
        Quaternion4 dg;
        dg.data[1]=-g.data[g.prd*1];
        dg.data[2]=-g.data[g.prd*2];
        dg.data[3]=-g.data[g.prd*3];
        dg.data[0]=g.data[g.prd*0];
        return dg;
    }

    /// Exponential map, using the inlined code from Eigen's conversion from axis/angle
    static Quaternion4 Expmap(const minivector& omega)
    {
        double theta2 = miniblas_vector_ddot(omega,omega);

        if (theta2 > std::numeric_limits<double>::epsilon())
        {
            double theta = sqrt(theta2);
            double ha = 0.5 * theta;
            double sinha_div_theta= (sin(ha) / theta);
            return Quaternion4(cos(ha), sinha_div_theta*omega.data[0],
                               sinha_div_theta*omega.data[1], sinha_div_theta*omega.data[2]);
        }
        else
        {
            // first order approximation sin(theta/2)/theta = 0.5
            return Quaternion4(1.0, 0.5*omega.data[0],
                               0.5*omega.data[1], 0.5*omega.data[2]);

        }
    }
    Quaternion4* ExpmapP(const minivector& omega)
    {
        double theta2 = miniblas_vector_ddot(omega,omega);

        if (theta2 > std::numeric_limits<double>::epsilon())
        {
            double theta = sqrt(theta2);
            double ha = 0.5 * theta;
            double sinha_div_theta= (sin(ha) / theta);
            return new Quaternion4(cos(ha), sinha_div_theta*omega.data[0],
                                   sinha_div_theta*omega.data[1], sinha_div_theta*omega.data[2]);
        }
        else
        {
            return new Quaternion4(1.0, 0.5*omega.data[0],
                                   0.5*omega.data[1], 0.5*omega.data[2]);

        }
    }

    static Quaternion4 Expmap(const minivector& omega,
                              minimatrix* H)
    {
        minimatrix_resize(H,3,3);
        minimatrix_memcpy(H,SO3::ExpmapDerivative(omega));
        double theta2 = miniblas_vector_ddot(omega,omega);

        if (theta2 > std::numeric_limits<double>::epsilon())
        {
            double theta = sqrt(theta2);
            double ha = 0.5 * theta;
            double sinha_div_theta= (sin(ha) / theta);
            // minivector vec = (sin(ha) / theta) * omega;
            return Quaternion4(cos(ha), sinha_div_theta*omega.data[0],
                               sinha_div_theta*omega.data[1], sinha_div_theta*omega.data[2]);
        }
        else
        {
            // first order approximation sin(theta/2)/theta = 0.5
            return Quaternion4(1.0, 0.5*omega.data[0],
                               0.5*omega.data[1], 0.5*omega.data[2]);

        }
    }

    /// We use our own Logmap, as there is a slight bug in Eigen
    static minivector Logmap(const Quaternion4& q)
    {

        // define these compile time constants to avoid std::abs:
        static const double twoPi = 2.0 * M_PI, NearlyOne = 1.0 - 1e-10,
                            NearlyNegativeOne = -1.0 + 1e-10;
        minivector omega(3);
        // const double qw = q.Q_.w();
        const double qw = q.data[0];
        double scale;
        // See Quaternion4-Logmap.nb in doc for Taylor expansions
        if (qw > NearlyOne)
        {
            // Taylor expansion of (angle / s) at 1
            // (2 + 2 * (1-qw) / 3) * q.vec();
            // omega = ( 8. / 3. - 2. / 3. * qw) * q.Q_.vec();
            scale=( 8. / 3. - 2. / 3. * qw) ;
        }
        else if (qw < NearlyNegativeOne)
        {
            // Taylor expansion of (angle / s) at -1
            // (-2 - 2 * (1 + qw) / 3) * q.vec();
            //  omega = (-8. / 3. - 2. / 3. * qw) * q.Q_.vec();
            scale=( -8. / 3. - 2. / 3. * qw) ;

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

            scale=angle / s ;
            //   omega = (angle / s) * q.Q_.vec();

        }
        omega.data[0] = scale * q.data[1];
        omega.data[1] = scale * q.data[2];
        omega.data[2] = scale * q.data[3];
        return omega;
    }

    static minivector Logmap(const Quaternion4& q, minimatrix* H)
    {

        // define these compile time constants to avoid std::abs:
        static const double twoPi = 2.0 * M_PI, NearlyOne = 1.0 - 1e-10,
                            NearlyNegativeOne = -1.0 + 1e-10;

        minivector omega(3);

        const double qw = q.data[0];
        double scale;
        // See Quaternion4-Logmap.nb in doc for Taylor expansions
        if (qw > NearlyOne)
        {
            // Taylor expansion of (angle / s) at 1
            // (2 + 2 * (1-qw) / 3) * q.vec();
            // omega = ( 8. / 3. - 2. / 3. * qw) * q.Q_.vec();
            scale=( 8. / 3. - 2. / 3. * qw) ;
        }
        else if (qw < NearlyNegativeOne)
        {
            // Taylor expansion of (angle / s) at -1
            // (-2 - 2 * (1 + qw) / 3) * q.vec();
            //  omega = (-8. / 3. - 2. / 3. * qw) * q.Q_.vec();
            scale=( -8. / 3. - 2. / 3. * qw) ;

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

            scale=angle / s ;
            //   omega = (angle / s) * q.Q_.vec();

        }
        omega.data[0] = scale * q.data[1];
        omega.data[1] = scale * q.data[2];
        omega.data[2] = scale * q.data[3];
        // *H = SO3::LogmapDerivative(omega);
        minimatrix mso3=SO3::LogmapDerivative(omega);
        minimatrix_memcpy(H,mso3);
        return omega;
    }

    /// @}
    /// @name Manifold traits
    /// @{

    static minivector Local(const Quaternion4& g, const Quaternion4& h)
    {
        Quaternion4 b = Between(g, h);
        minivector v = Logmap(b);
        return v;
    }

    static minivector Local(const Quaternion4& g, const Quaternion4& h,
                            minimatrix* H1, minimatrix* H2)
    {
        Quaternion4 b = Between(g, h, H1, H2);
        minimatrix D_v_b(3,3);
        minivector v = Logmap(b, &D_v_b);
        // *H1 = D_v_b * (*H1);
        minimatrix temp3(3,3);
        miniblas_dgemm(blasNoTrans,blasNoTrans,1.0,D_v_b,*H1,0.0,&temp3);
        minimatrix_memcpy(H1,temp3);

        // *H2 = D_v_b * (*H2);
        miniblas_dgemm(blasNoTrans,blasNoTrans,1.0,D_v_b,*H2,0.0,&temp3);
        minimatrix_memcpy(H2,temp3);
        return v;
    }
    virtual minimatrix* Retract(const minimatrix* mpose)
    {
        minivector v(mpose);
        return ExpmapP(v);

    }

    virtual minimatrix LocalCoordinates(const minimatrix* mpose,minimatrix* H1=NULL,minimatrix* H2=NULL) const
    {
        Quaternion4 v(minimatrix_get(mpose,0,0),minimatrix_get(mpose,1,0),minimatrix_get(mpose,2,0),minimatrix_get(mpose,3,0));
        return Logmap(v);
    }

    inline Quaternion4& operator=(const Quaternion4& rObj)
    {
        data[0]=rObj.data[0];
        data[1]=rObj.data[1];
        data[2]=rObj.data[2];
        data[3]=rObj.data[3];
        return *this;
    }

    ///@}

};
};

#endif // QUATERNION_H
