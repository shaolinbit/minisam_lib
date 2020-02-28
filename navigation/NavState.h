#ifndef NAVSTATE_H
#define NAVSTATE_H

/**
 * @file    NavState.h
 * @brief   Navigation state composing of attitude, position, and velocity
 **/

#pragma once

#include "../geometry/Pose3.h"
namespace minisam
{

/**
 * Navigation state: Pose (rotation, translation) + velocity
 */
class NavState:public minimatrix
{
public:


    typedef std::pair<minivector, minivector> PositionAndVelocity;

    /// @name Constructors
    /// @{

    /// Default constructor
    NavState() :
#ifdef USE_QUATERNIONS
        minimatrix(6,3)
#else
        minimatrix(5,3)
#endif
    {
        dimension=9;
        minimatrix_set_zero(this);
        data[0]=1.0;
#ifdef USE_QUATERNIONS

#else
        data[4]=1.0;
        data[8]=1.0;
#endif

    }
    ~NavState()
    {
    }

    /// Construct from attitude, position, velocity
    NavState(const Rot3& R, const minivector & t, const minivector & v) :
#ifdef USE_QUATERNIONS
        minimatrix(6,3)
#else
        minimatrix(5,3)
#endif
    {
        dimension=9;
#ifdef USE_QUATERNIONS
        minivector Rv=minivector_subcolumn(this,0,0,4);
        minivector_memcpy(&Rv,R);
        minivector t_=minimatrix_row(this,4);
        minivector_memcpy(&T_,t);
        minivector v_=minimatrix_row(this,5);
        minivector_memcpy(&v_,v);
#else
        minimatrix Rv=minimatrix_blockmatrix_var(this,0,0,3,3);
        minimatrix_memcpy(&Rv,R);
        minivector t_=minimatrix_row(this,3);
        minivector_memcpy(&t_,t);
        minivector v_=minimatrix_row(this,4);
        minivector_memcpy(&v_,v);
#endif
    }
    /// Construct from pose and velocity
    NavState(const Pose3& pose, const minivector & v):
#ifdef USE_QUATERNIONS
        minimatrix(6,3)
#else
        minimatrix(5,3)
#endif
    {
        dimension=9;
#ifdef USE_QUATERNIONS
        minivector Rv=minivector_subcolumn(this,0,0,4);
        minivector_memcpy(&pose_,pose.rotation());
        minivector t_=minimatrix_row(this,4);
        minivector_memcpy(&t_,pose.translation());
        minivector v_=minimatrix_row(this,5);
        minivector_memcpy(&v_,v);
#else
        minimatrix pose_=minimatrix_blockmatrix_var(this,0,0,4,3);
        minimatrix_memcpy(&pose_,pose);
        minivector v_=minimatrix_row(this,4);
        minivector_memcpy(&v_,v);
#endif
    }
    NavState(const minimatrix* pose, const minimatrix* v):
#ifdef USE_QUATERNIONS
        minimatrix(6,3)
#else
        minimatrix(5,3)
#endif
    {
        dimension=9;
#ifdef USE_QUATERNIONS
        minivector Rv=minivector_subcolumn(this,0,0,4);
        minivector subpose=minimatrix_subcolumn(pose,0,0,4);
        minivector_memcpy(&pose_,subpose);
        minivector t_=minimatrix_row(this,4);
        minivector t=minimatrix_subcolumn(pose,0,4,3);
        minivector_memcpy(&t_,t);
        minivector v_=minimatrix_row(this,5);
        minivector_memcpy(&v_,v);
#else
        minimatrix pose_=minimatrix_blockmatrix_var(this,0,0,4,3);
        minimatrix_memcpy(&pose_,pose);
        minivector v_=minimatrix_row(this,4);
        minimatrix_memcpy(&v_,v);
#endif
    }

    /// Construct from SO(3) and R^6
    NavState(const minimatrix& R, const minivector& tv) :
#ifdef USE_QUATERNIONS
        minimatrix(6,3)
#else
        minimatrix(5,3)
#endif
    {
        dimension=9;
#ifdef USE_QUATERNIONS
        minivector Rv=minivector_subcolumn(this,0,0,4);
        minivector_memcpy(&Rv,R);
        data[12]=minivector_get(&tv,0);
        data[13]=minivector_get(&tv,1);
        data[14]=minivector_get(&tv,2);
        data[15]=minivector_get(&tv,3);
        data[16]=minivector_get(&tv,4);
        data[17]=minivector_get(&tv,5);
#else
        minimatrix Rv=minimatrix_blockmatrix_var(this,0,0,3,3);
        minimatrix_memcpy(&Rv,R);
        data[9]=minivector_get(&tv,0);
        data[10]=minivector_get(&tv,1);
        data[11]=minivector_get(&tv,2);
        data[12]=minivector_get(&tv,3);
        data[13]=minivector_get(&tv,4);
        data[14]=minivector_get(&tv,5);
#endif
    }
    NavState(const minimatrix& state):
#ifdef USE_QUATERNIONS
        minimatrix(6,3)
#else
        minimatrix(5,3)
#endif
    {
        minimatrix_memcpy(this,state);
    }

    /// Named constructor with derivatives
    static NavState Create(const Rot3& R, const minivector& t, const minivector & v,
                           minimatrix* H1, minimatrix* H2,
                           minimatrix* H3);
    /// Named constructor with derivatives
    NavState FromPoseVelocity(const Pose3& pose, const minivector & vel,
                              minimatrix* H1, minimatrix* H2);

    /// @}
    /// @name Component Access
    /// @{


    Rot3 attitude(minimatrix*  H = NULL) const;
    //Rot3& attitude_();
    minivector  position(minimatrix* H = NULL) const;
    minivector  velocity(minimatrix* H = NULL) const;

    Pose3 pose() const
    {
        minimatrix pose_=minimatrix_blockmatrix(*this,0,0,4,3);
        return Pose3(&pose_);
    }

    /// @}
    /// @name Derived quantities
    /// @{

    /// Return rotation matrix. Induces computation in quaternion mode
    minimatrix R() const
    {
        return attitude().matrix();
    }
    /// Return quaternion. Induces computation in matrix mode
    Quaternion4 quaternion() const
    {
        return attitude().toQuaternion();
    }
    /// Return position as Vector3
    minivector  t() const
    {
#ifdef USE_QUATERNIONS
        minivector t_=minimatrix_row(*this,4);
#else
        minivector t_=minimatrix_row(*this,3);
#endif
        return t_;
    }
    /// Return velocity as Vector3. Computation-free.
    minivector v() const
    {
#ifdef USE_QUATERNIONS
        minivector v_=minimatrix_row(*this,5);
#else
        minivector v_=minimatrix_row(*this,4);
#endif
        return v_;
    }
    // Return velocity in body frame
    minivector  bodyVelocity(minimatrix* H = NULL) const;

    /// Return matrix group representation, in MATLAB notation:
    /// nTb = [nRb 0 n_t; 0 nRb n_v; 0 0 1]
    /// With this embedding in GL(3), matrix product agrees with compose
    minimatrix matrix() const;


    // Tangent space sugar.
    static  minivector dR(minivector& v)
    {
        return   minivector_subvector(v,0,3);
    }
    static  minivector dP(minivector& v)
    {
        return   minivector_subvector(v,3,3);
    }
    static  minivector dV(minivector& v)
    {
        return   minivector_subvector(v,6,3);
    }
    static  minivector dR(const minivector& v)
    {
        return   minivector_subvector(v,0,3);
    }
    static const minivector dP(const minivector& v)
    {
        return   minivector_subvector(v,3,3);
    }
    static  minivector dV(const minivector& v)
    {
        return   minivector_subvector(v,6,3);
    }

    /// retract with optional derivatives
    NavState retract(const minivector& v, //
                     minimatrix* H1 = NULL,minimatrix* H2 =
                         NULL)const ;
    virtual minimatrix* Retract(const minimatrix* mpose);
    /// localCoordinates with optional derivatives
    virtual minimatrix LocalCoordinates(const minimatrix* g, //
                                minimatrix* H1 = NULL, minimatrix* H2 =NULL) const;
    //virtual minimatrix LocalCoordinates(const minimatrix* mpose) const;

    /// @}
    /// @name Dynamics
    /// @{

    /// Integrate forward in time given angular velocity and acceleration in body frame
    /// Uses second order integration for position, returns derivatives except dt.
    NavState update(const minivector& b_acceleration, const minivector& b_omega,
                    const double dt, minimatrix* F, minimatrix* G1,
                    minimatrix* G2);

    /// Compute tangent space contribution due to Coriolis forces
    minivector coriolis(double dt, const minivector& omega, bool secondOrder = false,
                        minimatrix* H = NULL) const;

    /// Correct preintegrated tangent vector with our velocity and rotated gravity,
    /// taking into account Coriolis forces if omegaCoriolis is given.
    minivector correctPIM(const  minivector& pim, double dt, const minivector& n_gravity,
                          const minivector& omegaCoriolis, bool use2ndOrderCoriolisf =
                              false, minimatrix* H1 = NULL,
                          minimatrix* H2 =NULL) const;
    NavState& operator=(const NavState &rObj);

    void print() const;
};
};

#endif
