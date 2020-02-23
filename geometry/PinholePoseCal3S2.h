#ifndef PINHOLEPOSECAL3S2_H
#define PINHOLEPOSECAL3S2_H

/**
 * @file   PinholePose.h
 * @brief  Pinhole camera with known calibration
 * @author
 * @author
 * @date
 */
#pragma once

#include "../geometry/CalibratedCamera.h"
#include "../geometry/Cal3_S2.h"

namespace minisam
{

/**
 * A pinhole camera class that has a Pose3 and a *fixed* Calibration.
 * @addtogroup geometry
 * \nosubgrouping
 */
class PinholeBaseKCal3S2: public PinholeBase
{

//private:

    // Get dimensions of calibration type at compile time
    //static const int DimK = Cal3_S2::dimension;

public:

    /// @name Standard Constructors
    /// @{

    /** default constructor */
    PinholeBaseKCal3S2():PinholeBase()
    {

    }

    /** constructor with pose */
    explicit PinholeBaseKCal3S2(const Pose3& pose) :
        PinholeBase(pose)
    {
    }
    PinholeBaseKCal3S2(const minimatrix* m_memory)
    {
        size1=m_memory->size1;
    size2=m_memory->size2;
    prd=m_memory->prd;
    data=m_memory->data;
    owner=0;
    dimension=m_memory->dimension;
    }

    /// @}
    /// @name Advanced Constructors
    /// @{

    explicit PinholeBaseKCal3S2(const minivector &v) :
        PinholeBase(v)
    {
    }

    /// @}
    /// @name Standard Interface
    /// @{

    virtual ~PinholeBaseKCal3S2()
    {
    }

    /// return calibration
    virtual  Cal3_S2& calibration()
    {
        Cal3_S2 result;
        return result;
    }

    /// @}
    /// @name Transformations and measurement functions
    /// @{

    /// Project a point into the image and check depth
    std::pair<minivector, bool> projectSafe(const minivector& pw)
    {
        std::pair<minivector, bool> pn = PinholeBase::projectSafe(pw);
        // pn.first = calibration().uncalibrate(pn.first);
        minivector temp= calibration().uncalibrate(pn.first);
        minivector_memcpy(&pn.first,temp);
        return pn;
    }

    /** project a point from world coordinate to the image
     *  @param pw is a point in the world coordinates
     */
    minivector projectPoint(const minivector& pw)
    {
        minivector pn = PinholeBase::project2Point(pw); // project to normalized coordinates
        //return calibration().uncalibrate(pn); // uncalibrate to pixel coordinates
        minivector result=calibration().uncalibrate(pn);
        return result;
    }

    /** project a point from world coordinate to the image
     *  @param pw is a point at infinity in the world coordinates
     */
    minivector projectUnit(Unit3& pw)
    {
        const Unit3 pc = pose().rotation().unrotateUnit(pw); // convert to camera frame
        minivector pn = PinholeBase::ProjectUnit(pc); // project to normalized coordinates
        //return calibration().uncalibrate(pn);  // uncalibrate to pixel coordinates
        minivector result=calibration().uncalibrate(pn);  // uncalibrate to pixel coordinates
        return result;
    }


    /** Templated projection of a point (possibly at infinity) from world coordinate to the image
     *  @param pw is a 3D point or aUnit3 (point at infinity) in world coordinates
     *  @param Dpose is the Jacobian w.r.t. pose3
     *  @param Dpoint is the Jacobian w.r.t. point3
     *  @param Dcal is the Jacobian w.r.t. calibration
     */
    minivector _projectPoint(const minivector& pw,
                             minimatrix* Dpose,
                             minimatrix* Dpoint,
                             minimatrix* Dcal)
    {

        // project to normalized coordinates
        minivector pn = PinholeBase::project2Point(pw, Dpose, Dpoint);//

        // uncalibrate to pixel coordinates
        minimatrix Dpi_pn(2,2);
        minivector pi = calibration().uncalibrate(pn, Dcal,
                        Dpose || Dpoint ? &Dpi_pn : NULL);

        // If needed, apply chain rule
        if (Dpose!=NULL)
        {

            //   *Dpose = Dpi_pn * *Dpose;Dpose 2*6
            minimatrix temp1(2,6);
            minimatrix_memcpy(&temp1,*Dpose);
            miniblas_dgemm(blasNoTrans,blasNoTrans,1.0,Dpi_pn,temp1,0.0,Dpose);

        }
        if (Dpoint!=NULL)
        {
            // *Dpoint = Dpi_pn * *Dpoint;Dpoint 2*3
            minimatrix temp2(2,3);
            minimatrix_memcpy(&temp2,*Dpoint);
            miniblas_dgemm(blasNoTrans,blasNoTrans,1.0,Dpi_pn,*Dpoint,0.0,Dpoint);
        }




        return pi;
    }

    minivector _projectUnit(Unit3& pw,minimatrix* Dpose,
                            minimatrix* Dpoint,
                            minimatrix* Dcal)
    {

        // project to normalized coordinates
        minivector pn = PinholeBase::project2Unit(pw, Dpose, Dpoint);

        // uncalibrate to pixel coordinates
        minimatrix Dpi_pn(2,2);
        minivector pi = calibration().uncalibrate(pn, Dcal,
                        Dpose || Dpoint ? &Dpi_pn : NULL);

        // If needed, apply chain rule
        if (Dpose!=NULL)
        {
            //*Dpose = Dpi_pn * *Dpose;
            minimatrix temp1(2,6);
            minimatrix_memcpy(&temp1,*Dpose);
            miniblas_dgemm(blasNoTrans,blasNoTrans,1.0,Dpi_pn,temp1,0.0,Dpose);
        }
        if (Dpoint!=NULL)
        {
            //*Dpoint = Dpi_pn * *Dpoint;
            minimatrix temp2(2,3);
            minimatrix_memcpy(&temp2,*Dpoint);
            miniblas_dgemm(blasNoTrans,blasNoTrans,1.0,Dpi_pn,*Dpoint,0.0,Dpoint);
        }

        return pi;
    }

    /// project a 3D point from world coordinates into the image
    minivector projectPoint(const minivector& pw,minimatrix* Dpose,
                            minimatrix* Dpoint,
                            minimatrix* Dcal)
    {
        return _projectPoint(pw, Dpose, Dpoint, Dcal);
    }

    /// project a point at infinity from world coordinates into the image
    minivector projectUnit( Unit3& pw, minimatrix* Dpose,
                            minimatrix* Dpoint,
                            minimatrix* Dcal)
    {
        return _projectUnit(pw, Dpose, Dpoint, Dcal);
    }

    /// backproject a 2-dimensional point to a 3-dimensional point at given depth
    minivector backproject(const minivector& p, double depth)
    {
        minivector pn = calibration().calibrate(p);
        minivector result= pose().transform_from(backproject_from_camera(pn, depth));
        return result;
    }

    /// backproject a 2-dimensional point to a 3-dimensional point at infinity
    Unit3 backprojectPointAtInfinity(const minivector& p)
    {
        minivector pn = calibration().calibrate(p);
        Unit3 pc(pn.data[0], pn.data[1], 1.0); //by convention the last element is 1
        return pose().rotation().rotateUnit(pc);
    }

    /**
     * Calculate range to a landmark
     * @param point 3D location of landmark
     * @return range (double)
     */
    double range(const minivector& point,
                 minimatrix* Dcamera,
                 minimatrix* Dpoint) const
    {
        return pose().range(point, Dcamera, Dpoint);
    }

    /**
     * Calculate range to another pose
     * @param pose Other SO(3) pose
     * @return range (double)
     */
    double range(const Pose3& pose, minimatrix* Dcamera,
                 minimatrix*  Dpose) const
    {
        return this->pose().range(pose, Dcamera, Dpose);
    }

    /**
     * Calculate range to a CalibratedCamera
     * @param camera Other camera
     * @return range (double)
     */
    double range(const CalibratedCamera& camera, minimatrix*  Dcamera,
                 minimatrix*  Dother) const
    {
        return pose().range(camera.pose(), Dcamera, Dother);
    }

    /**
     * Calculate range to a PinholePoseK derived class
     * @param camera Other camera
     * @return range (double)
     */
    double range(const PinholeBaseKCal3S2& camera,
                 minimatrix* Dcamera,
                 minimatrix* Dother) const
    {
        return pose().range(camera.pose(), Dcamera, Dother);
    }

    ///@}


};
// end of class PinholeBaseK

/**
 * A pinhole camera class that has a Pose3 and a *fixed* Calibration.
 * Instead of using this class, one might consider calibrating the measurements
 * and using CalibratedCamera, which would then be faster.
 * @addtogroup geometry
 * \nosubgrouping
 */
class PinholePoseCal3S2:
#ifdef USE_QUATERNIONS
    public minivector//(12)
#else
    public minimatrix//(3,6)
#endif
{
public:

    /// @name Standard Constructors
    /// @{

    /** default constructor */
    PinholePoseCal3S2():
#ifdef USE_QUATERNIONS
        minivector(12)
#else
        minimatrix(3,6)
#endif
    {
        dimension=6;
#ifdef USE_QUATERNIONS
        data[0]=1.0;
        data[1]=0.0;
        data[2]=0.0;
        data[3]=0.0;
        data[4]=0.0;
        data[5]=0.0;
        data[6]=0.0;
        data[7]=1.0;
        data[8]=1.0;
        data[9]=0.0;
        data[10]=0.0;
        data[11]=0.0;

#else
        minimatrix_set_zero(this);
        data[0]=1.0;
        data[7]=1.0;
        data[14]=1.0;
        data[4]=1.0;
        data[10]=1.0;
#endif

    }

    /** constructor with pose, uses default calibration */
    explicit PinholePoseCal3S2(const Pose3& pose) ://PinholeBaseKCal3S2(pose), K_(new Cal3_S2())
#ifdef USE_QUATERNIONS
        minivector(12)
#else
        minimatrix(3,6)
#endif
    {
        dimension=6;
#ifdef USE_QUATERNIONS
        minivector qpose=minivector_subvector_var(this,0,7);
        minivector_memcpy(&qpose,pose);
        data[7]=1.0;
        data[8]=1.0;
        data[9]=0.0;
        data[10]=0.0;
        data[11]=0.0;
#else
        minimatrix qpose=minimatrix_blockmatrix_var(this,0,0,3,4);
        minimatrix_memcpy(&qpose,pose);
        data[4]=1.0;
        data[5]=0.0;
        data[10]=1.0;
        data[11]=0.0;
        data[16]=0.0;
        data[17]=0.0;
#endif

    }

    /** constructor with pose and calibration */
    PinholePoseCal3S2(const Pose3& pose,const Cal3_S2& K) :
#ifdef USE_QUATERNIONS
        minivector(12)
#else
        minimatrix(3,6)
#endif
    {
        dimension=6;
#ifdef USE_QUATERNIONS
        minivector qpose=minivector_subvector_var(this,0,7);
        minivector_memcpy(&qpose,pose);
        data[7]=K.data[0];
        data[8]=K.data[1];
        data[9]=K.data[2];
        data[10]=K.data[3];
        data[11]=K.data[4];
#else
        minimatrix qpose=minimatrix_blockmatrix_var(this,0,0,3,4);
        minimatrix_memcpy(&qpose,pose);
        data[4]=K.data[0];
        data[5]=K.data[3];
        data[10]=K.data[1];
        data[11]=K.data[4];
        data[16]=K.data[2];
        data[17]=0.0;
#endif
    }


    /// @}
    /// @name Named Constructors
    /// @{

    /**
     * Create a level camera at the given 2D pose and height
     * @param K the calibration
     * @param pose2 specifies the location and viewing direction
     * (theta 0 = looking in direction of positive X axis)
     * @param height camera height
     */
    static PinholePoseCal3S2 Level(const Cal3_S2& K,
                                   const Pose2& pose2, double height)
    {
        return PinholePoseCal3S2(PinholeBaseLevelPose(pose2, height), K);
    }

    /// PinholePose::level with default calibration
    static PinholePoseCal3S2 Level(const Pose2& pose2, double height)
    {
        return PinholePoseCal3S2::Level(Cal3_S2(), pose2, height);
    }

    /**
     * Create a camera at the given eye position looking at a target point in the scene
     * with the specified up direction vector.
     * @param eye specifies the camera position
     * @param target the point to look at
     * @param upVector specifies the camera up direction vector,
     *        doesn't need to be on the image plane nor orthogonal to the viewing axis
     * @param K optional calibration parameter
     */
    static PinholePoseCal3S2 Lookat(const minivector& eye, const minivector& target,
                                    const minivector& upVector,  const Cal3_S2& K =
                                        Cal3_S2())
    {
        return PinholePoseCal3S2(PinholeBaseLookatPose(eye, target, upVector), K);
    }

    /// @}
    /// @name Advanced Constructors
    /// @{

    /// Init from 6D vector
    explicit PinholePoseCal3S2(const minivector  &xi) ://PinholeBaseKCal3S2(v), K_(new Cal3_S2())
#ifdef USE_QUATERNIONS
        minivector(12)
#else
        minimatrix(3,6)
#endif
    {
        dimension=6;
        minivector omega=minivector_subvector(xi,0,3);
        minivector v=minivector_subvector(xi,3,3);
        Rot3 R = Rot3::Expmap(omega);
        double theta2 = miniblas_vector_ddot(omega,omega);

        Pose3* result;

        if (theta2 > std::numeric_limits<double>::epsilon())
        {
            double dotomegav=miniblas_vector_ddot(omega,v);
            minivector t_parallel(3);
            minivector_scale_vec(&t_parallel,dotomegav,omega);
            minivector omega_cross_v =cross3d(omega,v);

            minivector t=R.multiplyvector(omega_cross_v);
            minivector_scale(&t,-1.0);
            minivector_add(&t,omega_cross_v);
            minivector_add(&t,t_parallel);
            minivector_scale(&t,1/theta2);
            result=new Pose3(R,t);
        }
        else
        {
            result=new Pose3(R, v);
        }

#ifdef USE_QUATERNIONS
        minivector qpose=minivector_subvector_var(this,0,7);
        minivector_memcpy(&qpose,result);
        data[7]=1.0;
        data[8]=1.0;
        data[9]=0.0;
        data[10]=0.0;
        data[11]=0.0;
#else
        minimatrix qpose=minimatrix_blockmatrix_var(this,0,0,3,4);
        minimatrix_memcpy(&qpose,result);
        data[4]=1.0;
        data[5]=0.0;
        data[10]=1.0;
        data[11]=0.0;
        data[16]=0.0;
        data[17]=0.0;
#endif

    }

    /// Init from Vector and calibration
    PinholePoseCal3S2(const minivector& xi, const minivector &K) ://PinholeBaseKCal3S2(v), K_(new Cal3_S2(K))
#ifdef USE_QUATERNIONS
        minivector(12)
#else
        minimatrix(3,6)
#endif
    {
        dimension=6;
        minivector omega=minivector_subvector(xi,0,3);
        minivector v=minivector_subvector(xi,3,3);
        Rot3 R = Rot3::Expmap(omega);
        double theta2 = miniblas_vector_ddot(omega,omega);

        Pose3* result;

        if (theta2 > std::numeric_limits<double>::epsilon())
        {
            double dotomegav=miniblas_vector_ddot(omega,v);
            minivector t_parallel(3);
            minivector_scale_vec(&t_parallel,dotomegav,omega);
            minivector omega_cross_v =cross3d(omega,v);

            minivector t=R.multiplyvector(omega_cross_v);
            minivector_scale(&t,-1.0);
            minivector_add(&t,omega_cross_v);
            minivector_add(&t,t_parallel);
            minivector_scale(&t,1/theta2);
            result=new Pose3(R,t);
        }
        else
        {
            result=new Pose3(R, v);
        }

#ifdef USE_QUATERNIONS
        minivector qpose=minivector_subvector_var(this,0,7);
        minivector_memcpy(&qpose,result);
        minivector qk=minivector_subvector_var(this,7,5);
        minivector_memcpy(&qk,K);
#else
        minimatrix qpose=minimatrix_blockmatrix_var(this,0,0,3,4);
        minimatrix_memcpy(&qpose,result);
        data[4]=minivector_get(&K,0);
        data[5]=minivector_get(&K,3);
        data[10]=minivector_get(&K,1);
        data[11]=minivector_get(&K,4);
        data[16]=minivector_get(&K,2);
        data[17]=0.0;
#endif
    }

    /*
    /// stream operator
    friend std::ostream& operator<<(std::ostream &os, const PinholePoseCal3S2& camera)
    {
        os << "{R: " << camera.pose().rotation().rpy().transpose();
        os << ", t: " << camera.pose().translation().transpose();
        if (!camera.K_)
            os << ", K: none";
        else
            os << ", K: " << *camera.K_;
        os << "}";
        return os;
    }*/
    /// @}
    /// @name Standard Interface
    /// @{

    ~PinholePoseCal3S2()
    {
    }

    /// return  pointer to calibration
    Cal3_S2 calibration() const
    {
        //return K_;
#ifdef USE_QUATERNIONS
        minivector qk=minivector_subvector(*this,7,5);
        return Cal3_S2(qk);
#else
        return Cal3_S2(data[4],data[10],data[16],data[5],data[11]);
#endif
    }

    /// return calibration
    //virtual  Cal3_S2& calibration()
    //{
    //    return *K_;
    // }

    /** project a point from world coordinate to the image, 2 derivatives only
     *  @param pw is a point in world coordinates
     *  @param Dpose is the Jacobian w.r.t. the whole camera (really only the pose)
     *  @param Dpoint is the Jacobian w.r.t. point3
     */
    minivector project2Point(const minivector& pw, minimatrix* Dpose,
                             minimatrix* Dpoint)
    {
#ifdef USE_QUATERNIONS
        minivector qpose=minivector_subvector(*this,0,7);
#else
        minimatrix qpose=minimatrix_blockmatrix(*this,0,0,3,4);
#endif // USE_QUATERNIONS
        PinholeBaseKCal3S2 result(&qpose);
        //minimatrix_memcpy(&result,qpose);
        return result.projectPoint(pw, Dpose, Dpoint,NULL);
    }

    /// project2 version for point at infinity
    minivector project2Unit(Unit3& pw, minimatrix* Dpose,
                            minimatrix* Dpoint)
    {
#ifdef USE_QUATERNIONS
        minivector qpose=minivector_subvector(*this,0,7);
#else
        minimatrix qpose=minimatrix_blockmatrix(*this,0,0,3,4);
#endif // USE_QUATERNIONS
        PinholeBaseKCal3S2 result(&qpose);
        return result.projectUnit(pw, Dpose, Dpoint,NULL);
    }

    /// @}
    /// @name Manifold
    /// @{

    int dim() const
    {
        return 6;
    }

    /// move a cameras according to d
    virtual minimatrix* Retract(const minimatrix* mpose)
    {
#ifdef USE_QUATERNIONS
        minivector qk=minivector_subvector(*this,7,5);
        return new PinholePoseCal3S2(Pose3::ChartAtOrigin::retract(*mpose), Cal3_S2(&qk));
#else
        return new PinholePoseCal3S2(Pose3::ChartAtOrigin::retract(*mpose),  Cal3_S2(data[4],data[10],data[16],data[5],data[11]));
#endif
    }
    virtual minimatrix LocalCoordinates(const minimatrix* mpose) const
    {

#ifdef USE_QUATERNIONS
        minivector qpose=minivector_subvector(*this,0,7);
        Pose3 result(qpose);
        return Pose3::ChartAtOrigin::Local(result);
#else
        minimatrix qpose=minimatrix_blockmatrix(*this,0,0,3,4);
        Pose3 result;
        minimatrix_memcpy(&result,qpose);
        return Pose3::ChartAtOrigin::Local(result);
#endif // USE_QUATERNIONS

    }
    /*
    PinholePoseCal3S2 retract(const minivector& d) const
    {
        return PinholePoseCal3S2(Pose3::ChartAtOrigin::retract(d), K_);
    }


    /// return canonical coordinate
    minivector LocalCoordinates(const PinholePoseCal3S2& p) const
    {
        return Pose3::ChartAtOrigin::Local(p.PinholeBaseKCal3S2::pose());
    }
    */

    /// for Canonical
    static PinholePoseCal3S2 identity()
    {
        return PinholePoseCal3S2(); // assumes that the default constructor is valid
    }

    /// @}
};

};
// end of class PinholePose
#endif // PINHOLEPOSE_H
