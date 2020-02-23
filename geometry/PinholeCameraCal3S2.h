#ifndef PINHOLECAMERACAL3S2_H
#define PINHOLECAMERACAL3S2_H


/**
 * @file PinholeCamera.h
 * @brief Base class for all pinhole cameras
 * @author
 * @date
 */

#pragma once

#include "../geometry/PinholePoseCal3S2.h"

namespace minisam
{

/**
 * A pinhole camera class that has a Pose3 and a Calibration.
 * Use PinholePose if you will not be optimizing for Calibration
 * @addtogroup geometry
 * \nosubgrouping
 */
class PinholeCameraCal3S2:// public PinholeBaseKCal3S2
#ifdef USE_QUATERNIONS
    public minivector//(12)
#else
    public minimatrix//(6,3)
#endif
{
//private:

    //  Cal3_S2 K_; ///< Calibration, part of class now

    // Get dimensions of calibration type at compile time
    //  static const int DimK = Cal3_S2::dimension;

public:

    // enum
    // {
    //      dimension = 6 + DimK
//   }; ///< Dimension depends on calibration

    /// @name Standard Constructors
    /// @{

    /** default constructor */
    PinholeCameraCal3S2():
#ifdef USE_QUATERNIONS
        minivector(12)
#else
        minimatrix(6,3)
#endif
    {
        dimension=11;
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
        data[4]=1.0;
        data[8]=1.0;
        data[12]=1.0;
        data[13]=1.0;
#endif

    }
    /** constructor with pose */
    explicit PinholeCameraCal3S2(const Pose3& pose) :
#ifdef USE_QUATERNIONS
        minivector(12)
#else
        minimatrix(6,3)
#endif
    {
        dimension=11;
#ifdef USE_QUATERNIONS
        minivector qpose=minivector_subvector_var(this,0,7);
        minivector_memcpy(&qpose,pose);
        data[7]=1.0;
        data[8]=1.0;
        data[9]=0.0;
        data[10]=0.0;
        data[11]=0.0;
#else
        minimatrix qpose=minimatrix_blockmatrix_var(this,0,0,4,3);
        minimatrix_memcpy(&qpose,pose);
        data[12]=1.0;
        data[13]=1.0;
        data[14]=0.0;
        data[15]=0.0;
        data[16]=0.0;
        data[17]=0.0;
#endif

    }


    /** constructor with pose and calibration */
    PinholeCameraCal3S2(const Pose3& pose, const Cal3_S2& K) :
#ifdef USE_QUATERNIONS
        minivector(12)
#else
        minimatrix(6,3)
#endif    //PinholeBaseKCal3S2(pose), K_(K)
    {
        dimension=11;
#ifdef USE_QUATERNIONS
        minivector qpose=minivector_subvector_var(this,0,7);
        minivector_memcpy(&qpose,pose);
        data[7]=K.data[0];
        data[8]=K.data[1];
        data[9]=K.data[2];
        data[10]=K.data[3];
        data[11]=K.data[4];
#else
        minimatrix qpose=minimatrix_blockmatrix_var(this,0,0,4,3);
        minimatrix_memcpy(&qpose,pose);
        data[12]=K.data[0];
        data[13]=K.data[1];
        data[14]=K.data[2];
        data[15]=K.data[3];
        data[16]=K.data[4];
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
    static PinholeCameraCal3S2 Level(const Cal3_S2 &K, const Pose2& pose2,
                                     double height)
    {
        return PinholeCameraCal3S2(PinholeBaseLevelPose(pose2, height), K);
    }

    /// PinholeCamera::level with default calibration
    static PinholeCameraCal3S2 Level(const Pose2& pose2, double height)
    {
        return PinholeCameraCal3S2::Level(Cal3_S2(), pose2, height);
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
    static PinholeCameraCal3S2 Lookat(const minivector& eye, const minivector& target,
                                      const minivector& upVector, const Cal3_S2& K = Cal3_S2())
    {
        return PinholeCameraCal3S2(PinholeBaseLookatPose(eye, target, upVector), K);
    }

    // Create PinholeCamera, with derivatives
    static PinholeCameraCal3S2 Create(const Pose3& pose, const Cal3_S2 &K,
                                      minimatrix* H1, //
                                      minimatrix* H2)
    {
        if (H1!=NULL)//11*6
        {
            //*H1<<minimatrix::Identity(6,6),minimatrix::Zero(DimK,6);
            minimatrix block1=minimatrix_blockmatrix_var(H1,0,0,6,6);
            minimatrix_set_identity(&block1);
            minimatrix block2=minimatrix_blockmatrix_var(H1,6,0,5,6);
            minimatrix_set_zero(&block2);
        }
        if (H2!=NULL)
        {
            //*H2 << minimatrix::Zero(6,DimK),minimatrix::Identity(DimK,DimK);
            minimatrix block3=minimatrix_blockmatrix_var(H2,0,0,6,5);
            minimatrix_set_zero(&block3);
            minimatrix block4=minimatrix_blockmatrix_var(H2,6,0,5,5);
            minimatrix_set_identity(&block4);

        }
        return PinholeCameraCal3S2(pose,K);
    }

    /// @}
    /// @name Advanced Constructors
    /// @{

    /// Init from vector, can be 6D (default calibration) or dim
    explicit PinholeCameraCal3S2(const minivector&xi) :
#ifdef USE_QUATERNIONS
        minivector(12)
#else
        minimatrix(6,3)
#endif
    {
        dimension=11;
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
        if(v.size1==11)
        {
            data[7]=minivector_get(&v,7);
            data[8]=minivector_get(&v,8);
            data[9]=minivector_get(&v,9);
            data[10]=minivector_get(&v,10);
            data[11]=minivector_get(&v,11);
        }
#else
        minimatrix qpose=minimatrix_blockmatrix_var(this,0,0,4,3);
        minimatrix_memcpy(&qpose,result);
        if(v.size1==11)
        {
            data[12]=minivector_get(&v,7);
            data[13]=minivector_get(&v,8);
            data[14]=minivector_get(&v,9);
            data[15]=minivector_get(&v,10);
            data[16]=minivector_get(&v,11);
            data[17]=0.0;
        }
#endif

    }
    /// Init from Vector and calibration
    PinholeCameraCal3S2(const minivector& xi, const minivector&K) ://PinholeBaseKCal3S2(v), K_(K)
#ifdef USE_QUATERNIONS
        minivector(12)
#else
        minimatrix(6,3)
#endif
    {
        dimension=11;
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
        minimatrix qpose=minimatrix_blockmatrix_var(this,0,0,4,3);
        minimatrix_memcpy(&qpose,result);
        data[12]=minivector_get(&K,0);
        data[13]=minivector_get(&K,1);
        data[14]=minivector_get(&K,2);
        data[15]=minivector_get(&K,3);
        data[16]=minivector_get(&K,4);
        data[17]=0.0;
#endif
    }
    /// @}
    /// @name Standard Interface
    /// @{

    ~PinholeCameraCal3S2()
    {
    }

    /// return pose
    Pose3 pose() const
    {
#ifdef USE_QUATERNIONS
        minivector f=minivector_subvector(*this,0,7);
        return Pose3(&f);
#else
        minimatrix f=minimatrix_blockmatrix(*this,0,0,4,3);
        return Pose3(&f);
#endif

    }

    /// return pose, with derivative
     Pose3 getPose(minimatrix* H) const
    {
        if (H!=NULL)
        {
            minimatrix_set_zero(H);
            minimatrix block=minimatrix_blockmatrix_var(H,0,0,6,6);
            minimatrix_set_identity(&block);
        }
        return pose();
    }

    /// return calibration
     Cal3_S2 calibration()
    {
#ifdef USE_QUATERNIONS
        minivector qk=minivector_subvector(*this,7,5);
        return Cal3_S2(&qk);
#else
        return Cal3_S2(data[12],data[13],data[14],data[15],data[16]);
#endif
    }

    /// @}
    /// @name Manifold
    /// @{

    /// @deprecated
    int dim() const
    {
        return dimension;
    }



    /// move a cameras according to d
    virtual minimatrix* Retract(const minimatrix* mpose)
    {
        minivector mm(mpose);
        Cal3_S2 cali=calibration();
        minivector rsub=minivector_subvector(mm,mm.size1-cali.dim(),cali.dim());
        minimatrix* bb=cali.Retract(&rsub);
        Cal3_S2 bres(minimatrix_get(bb,0,0),minimatrix_get(bb,1,0),minimatrix_get(bb,2,0),
                     minimatrix_get(bb,3,0),minimatrix_get(bb,4,0));
        return new PinholeCameraCal3S2(Pose3::ChartAtOrigin::retract(minivector_subvector(mm,0,6)),bres);
    }

    virtual minimatrix LocalCoordinates(const minimatrix* T2) const
    {
        minivector d(dimension);
        minivector dhead=minivector_subvector_var(&d,0,6);

        PinholeCameraCal3S2 T2c;
        minimatrix_memcpy(&T2c,*T2);

        minivector head=Pose3::ChartAtOrigin::Local(T2c.pose());
        minivector_memcpy(&dhead,head);
        minivector dtail=minivector_subvector_var(&d,6,5);

        Cal3_S2 bb=T2c.calibration();
#ifdef USE_QUATERNIONS
        minivector qk=minivector_subvector(*this,7,5);
        Cal3_S2 cali(qk);
#else
        Cal3_S2 cali(data[12],data[13],data[14],data[15],data[16]);
#endif

        minimatrix lresult=cali.LocalCoordinates(&bb);
        minimatrix_memcpy(&dtail,lresult);
        return d;
    }
    /// for Canonical
    static PinholeCameraCal3S2 identity()
    {
        return PinholeCameraCal3S2(); // assumes that the default constructor is valid
    }

    /// @}
    /// @name Transformations and measurement functions
    /// @{

    /** Templated projection of a 3D point or a point at infinity into the image
     *  @param pw either a Point3 or a Unit3, in world coordinates
     */
    minivector _project2Point(const minivector& pw, minimatrix*  Dcamera,
                              minimatrix* Dpoint)
    {
        minimatrix Dpose(2,6);
        minimatrix Dcal(2,5);
#ifdef USE_QUATERNIONS
        minivector qpose=minivector_subvector(*this,0,7);
#else
        minimatrix qpose=minimatrix_blockmatrix(*this,0,0,4,3);
#endif // USE_QUATERNIONS

        PinholeBaseKCal3S2 rb(&qpose);
        minivector pi = rb.projectPoint(pw, Dcamera ? &Dpose : NULL, Dpoint,
                                        Dcamera ? &Dcal : NULL);
        if (Dcamera!=NULL)
        {
            minimatrix dblock1=minimatrix_blockmatrix_var(Dcamera,0,0,2,6);
            minimatrix_memcpy(&dblock1,Dpose);
            minimatrix dblock2=minimatrix_blockmatrix_var(Dcamera,6,0,2,5);
            minimatrix_memcpy(&dblock2,Dcal);

        }

        return pi;
    }
    minivector _project2Unit( Unit3& pw, minimatrix*  Dcamera,
                              minimatrix* Dpoint)
    {
        minimatrix Dpose(2,6);
        minimatrix Dcal(2,5);

#ifdef USE_QUATERNIONS
        minivector qpose=minivector_subvector(*this,0,7);
#else
        minimatrix qpose=minimatrix_blockmatrix(*this,0,0,4,3);
#endif // USE_QUATERNIONS

        PinholeBaseKCal3S2 rb(&qpose);

        minivector pi = rb.projectUnit(pw, Dcamera ? &Dpose : NULL, Dpoint,
                                       Dcamera ? &Dcal : NULL);
        if (Dcamera!=NULL)
        {
            minimatrix_resize(Dcamera,2,11);
            minimatrix dblock1=minimatrix_blockmatrix_var(Dcamera,0,0,2,6);
            minimatrix_memcpy(&dblock1,Dpose);
            minimatrix dblock2=minimatrix_blockmatrix_var(Dcamera,6,0,2,5);
            minimatrix_memcpy(&dblock2,Dcal);

        }
        return pi;
    }

    /// project a 3D point from world coordinates into the image
    minivector project2Point(const minivector& pw,  minimatrix*  Dcamera,
                             minimatrix*  Dpoint)
    {
        return _project2Point(pw, Dcamera, Dpoint);
    }

    /// project a point at infinity from world coordinates into the image
    minivector project2Unit( Unit3& pw, minimatrix*  Dcamera,
                             minimatrix*  Dpoint)
    {
        return _project2Unit(pw, Dcamera, Dpoint);
    }

       /// project a 3D point from world coordinates into the image
    minivector projectPoint(const minivector& pw,minimatrix* Dpose,
                            minimatrix* Dpoint,
                            minimatrix* Dcal)
    {

    minimatrix Rt; // calculated by transform_to if needed
     minivector q = pose().transform_toDPoint(pw, &Rt);
#ifdef THROW_CHEIRALITY_EXCEPTION
    if (q.z() <= 0)
        throw CheiralityException();
#endif

     double d = 1.0 / q.data[2];
     const double u = q.data[0] * d, v =q.data[1] * d;
    minivector pn =  minivector_2dim(u, v);


    if (Dpose!=NULL || Dpoint!=NULL) {
    const double d = 1.0 / q.data[2];
      if (Dpose!=NULL)
    {
    minimatrix_resize(Dpose,2,6);
    minimatrix_memcpy(Dpose,PinholeBaseDpose(pn, d));
    }
      if (Dpoint!=NULL)
   {
   minimatrix_resize(Dpoint,2,3);
   minimatrix_memcpy(Dpoint,PinholeBaseDpoint(pn, d, Rt));
   if(DEBUGSTATE)
   {
   std::cout<<"Dpoint"<<std::endl;
   minimatrix_print(Dpoint);
   std::cout<<std::endl;
   }
   }
 }

        minimatrix Dpi_pn(2,2);
        minivector pi = calibration().uncalibrate(pn, Dcal,
                        Dpose || Dpoint ? &Dpi_pn : NULL);

        if (Dpose!=NULL)
        {
            minimatrix temp1(2,6);
            minimatrix_memcpy(&temp1,*Dpose);
            miniblas_dgemm(blasNoTrans,blasNoTrans,1.0,Dpi_pn,temp1,0.0,Dpose);

        }
        if (Dpoint!=NULL)
        {
            minimatrix temp2(2,3);
            minimatrix_memcpy(&temp2,Dpoint);
            miniblas_dgemm(blasNoTrans,blasNoTrans,1.0,Dpi_pn,temp2,0.0,Dpoint);
            if(DEBUGSTATE)
            {
             std::cout<<"Dpoint"<<std::endl;
            minimatrix_print(Dpoint);
            std::cout<<std::endl;
            }
        }
        return pi;
    }

    /**
     * Calculate range to a landmark
     * @param point 3D location of landmark
     * @return range (double)
     */
    double range(const minivector& point, minimatrix* Dcamera,
                 minimatrix* Dpoint) const
    {
        minimatrix Dpose_(1,6);

        double result = this->pose().range(point, Dcamera ? &Dpose_ : NULL, Dpoint);
        if (Dcamera!=NULL)
        {
            minimatrix db1=minimatrix_blockmatrix_var(Dcamera,0,0,1,6);
            minimatrix_memcpy(&db1,Dpose_);
            minimatrix db2=minimatrix_blockmatrix_var(Dcamera,6,0,1,5);
            minimatrix_set_zero(&db2);

        }
        return result;
    }

    /**
     * Calculate range to another pose
     * @param pose Other SO(3) pose
     * @return range (double)
     */
    double range(const Pose3& pose, minimatrix* Dcamera,
                 minimatrix* Dpose) const
    {
        minimatrix Dpose_(1,6);

        double result = this->pose().range(pose, Dcamera ? &Dpose_ : NULL, Dpose);
        if (Dcamera!=NULL)
        {
            minimatrix db1=minimatrix_blockmatrix_var(Dcamera,0,0,1,6);
            minimatrix_memcpy(&db1,Dpose_);
            minimatrix db2=minimatrix_blockmatrix_var(Dcamera,6,0,1,5);
            minimatrix_set_zero(&db2);

        }
        return result;
    }

    /**
     * Calculate range to another camera
     * @param camera Other camera
     * @return range (double)
     */
    double range(const PinholeCameraCal3S2& camera,
                 minimatrix*  Dcamera,
                 minimatrix*  Dother) const
    {
        minimatrix Dcamera_(1,6), Dother_(1,6);

        double result = this->pose().range(camera.pose(), Dcamera ? &Dcamera_ : NULL,
                                           Dother ? &Dother_ : NULL);
        if (Dcamera!=NULL)
        {
            minimatrix db1=minimatrix_blockmatrix_var(Dcamera,0,0,1,6);
            minimatrix_memcpy(&db1,Dcamera_);
            minimatrix db2=minimatrix_blockmatrix_var(Dcamera,6,0,1,5);
            minimatrix_set_zero(&db2);

        }
        if (Dother!=NULL)
        {
            minimatrix_set_zero(Dother);
            minimatrix db3=minimatrix_blockmatrix_var(Dother,0,0,1,6);
            minimatrix_memcpy(&db3,Dother_);
        }
        return result;
    }

    /**
     * Calculate range to a calibrated camera
     * @param camera Other camera
     * @return range (double)
     */
    double range(const CalibratedCamera& camera,
                 minimatrix*  Dcamera,
                 minimatrix*   Dother) const
    {
        return range(camera.pose(), Dcamera, Dother);
    }

    ///@}

};
};
#endif // PINHOLECAMERA_H
