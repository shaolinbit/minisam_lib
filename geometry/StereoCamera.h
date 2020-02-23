/**
 * @file    StereoCamera.h
 * @brief   A Rectified Stereo Camera
 */

#pragma once

#include "../geometry/Cal3_S2Stereo.h"
#include "../geometry/Pose3.h"
#include "../geometry/StereoPoint2.h"

namespace minisam
{


/**
 * A stereo camera class, parameterize by left camera pose and stereo calibration
 * @addtogroup geometry
 */
class  StereoCamera:
#ifdef USE_QUATERNIONS
    public minivector//(13)
#else
    public minimatrix//(3,6)
#endif
{

public:

    /**
     *  Some classes template on either PinholeCamera or StereoCamera,
     *  and this typedef informs those classes what "project" returns.
     */

public:

    /// @name Standard Constructors
    /// @{

    /// Default constructor allocates a calibration!
    StereoCamera() :
#ifdef USE_QUATERNIONS
        minivector(13)
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
        data[12]=0.0;

#else
        minimatrix_set_zero(this);
        data[0]=1.0;
        data[7]=1.0;
        data[14]=1.0;
        data[4]=1.0;
        data[10]=1.0;
#endif

    }

    StereoCamera(const minimatrix* mcopy) :
#ifdef USE_QUATERNIONS
        minivector(13)
#else
        minimatrix(3,6)
#endif
    {
        dimension=6;
#ifdef USE_QUATERNIONS
        if(mcopy->size1!=13||mcopy->size2!=1)
        {
            throw std::invalid_argument("size not right");
        }
#else
        if(mcopy->size1!=3||mcopy->size2!=6)
        {
            throw std::invalid_argument("size not right");
        }
#endif
        prd=mcopy->prd;
        data=mcopy->data;
        owner=0;
    }

    /// Construct from pose and shared calibration
    StereoCamera(const Pose3& leftCamPose, const Cal3_S2Stereo& K);


    /// Return shared pointer to calibration
    Cal3_S2Stereo calibration() const
    {
#ifdef USE_QUATERNIONS
        minivector qk=minivector_subvector(*this,7,6);
        return Cal3_S2Stereo(&qk);
#else
        return Cal3_S2Stereo(data[4],data[10],data[16],data[5],data[11],data[17]);
#endif

    }

    /// @}


    /// @name Manifold
    /// @{

    /// Dimensionality of the tangent space
    inline size_t dim() const
    {
        return 6;
    }

    virtual minimatrix* Retract(const minimatrix* mpose)
    {
        Pose3 pp(pose());
        minimatrix* bp=pp.Retract(mpose);
        return new StereoCamera(Pose3(bp),calibration());
    }

    virtual minimatrix LocalCoordinates(const minimatrix* mpose) const
    {
        StereoCamera t2(mpose);
        Pose3 t2pose=t2.pose();
        return pose().LocalCoordinates(&t2pose);

    }


    /// @}
    /// @name Standard Interface
    /// @{

    /// pose
    Pose3 pose() const
    {
#ifdef USE_QUATERNIONS
        minivector f=minivector_subvector(*this,0,4);
        return Pose3(f);
#else
        minimatrix f=minimatrix_blockmatrix(*this,0,0,3,4);
        return Pose3(&f);
#endif
    }

    /// baseline
    double baseline() const
    {
#ifdef USE_QUATERNIONS
        return data[12];
#else
        return data[17];
#endif
    }

    /// Project 3D point to StereoPoint2 (uL,uR,v)
    StereoPoint2 project(const minivector& point) const;

    /** Project 3D point and compute optional derivatives
     * @param H1 derivative with respect to pose
     * @param H2 derivative with respect to point
     */
    StereoPoint2 project2(const minivector& point, minimatrix* H1=NULL,//3*6,
                          minimatrix* H2=NULL//3*3
                         ) const;

    /// back-project a measurement
    minivector backproject(const StereoPoint2& z) const;

    /** Back-project the 2D point and compute optional derivatives
     * @param H1 derivative with respect to pose
     * @param H2 derivative with respect to point
     */
    minivector backproject2(const StereoPoint2& z,
                            minimatrix* H1=NULL,//3*6,
                            minimatrix* H2=NULL//3*3
                           ) const;

    /// @}


};

};
