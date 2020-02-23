/*
 * @file EssentialMatrix.h
 * @brief EssentialMatrix class
 */

#pragma once

#include "../geometry/Pose3.h"
#include "..//geometry/Unit3.h"
#include "../mat/MatCal.h"
#include "../mat/Matrix.h"

#include <iosfwd>
#include <string>

namespace minisam
{

/**
 * An essential matrix is like a Pose3, except with translation up to scale
 * It is named after the 3*3 matrix aEb = [aTb]x aRb from computer vision,
 * but here we choose instead to parameterize it as a (Rot3,Unit3) pair.
 * We can then non-linearly optimize immediately on this 5-dimensional manifold.
 */
class  EssentialMatrix:public minimatrix
{
public:
    /// Static function to convert Point2 to homogeneous coordinates
    static minivector Homogeneous(const minivector& p)
    {
        minivector result(3);
        result.data[0]=minivector_get(&p,0);
        result.data[1]=minivector_get(&p,1);
        result.data[2]=1.0;
        return result;
    }

    /// @name Constructors and named constructors
    /// @{

    /// Default constructor
    EssentialMatrix() ://E_(t_.skew())
#ifdef USE_QUATERNIONS
        minimatrix(5,4)
#else
        minimatrix(7,3)
#endif
    {
        dimension=5;
#ifdef USE_QUATERNIONS
        data[0]=0.0;
        data[1]=0.0;
        data[2]=0.0;
        data[3]=0.0;

        data[4]=0.0;
        data[5]=0.0;
        data[6]=-1.0;
        data[7]=0.0;

        data[8]=0.0;
        data[9]=1.0;
        data[10]=0.0;
        data[11]=0.0;

        data[12]=1.0;
        data[13]=0.0;
        data[14]=0.0;

        data[15]=0.0;

        data[16]=1.0;
        data[17]=0.0;
        data[18]=0.0;
        data[19]=0.0;
#else
        data[0]=0.0;
        data[1]=0.0;
        data[2]=0.0;

        data[3]=0.0;
        data[4]=0.0;
        data[5]=-1.0;

        data[6]=0.0;
        data[7]=1.0;
        data[8]=0.0;

        data[9]=1.0;
        data[10]=0.0;
        data[11]=0.0;

        data[12]=1.0;
        data[13]=0.0;
        data[14]=0.0;

        data[15]=0.0;
        data[16]=1.0;
        data[17]=0.0;

        data[18]=0.0;
        data[19]=0.0;
        data[20]=1.0;

#endif
    }

    /// Construct from rotation and translation
    EssentialMatrix(const Rot3& aRb, const Unit3& aTb) :
#ifdef USE_QUATERNIONS
        minimatrix(5,4)
#else
        minimatrix(7,3)
#endif
    {
#ifdef USE_QUATERNIONS
        minimatrix mblock=minimatrix_blockmatrix_var(this,0,0,3,3);
        miniblas_dgemm(blasNoTrans,blasNoTrans,1.0,aTb.skew(),aRb.matrix(),0.0,&mblock);

        data[3]=0.0;
        data[7]=0.0;
        data[11]=0.0;

        data[12]=aTb.data[0];
        data[13]=aTb.data[1];
        data[14]=aTb.data[2];
        data[15]=0.0;

        minivector subq=minimatrix_row(this,4);
        minivector_memcpy(&subq,aRb);
#else
        minimatrix mblock=minimatrix_blockmatrix_var(this,0,0,3,3);
        miniblas_dgemm(blasNoTrans,blasNoTrans,1.0,aTb.skew(),aRb.matrix(),0.0,&mblock);

        data[9]=aTb.data[0];
        data[10]=aTb.data[1];
        data[11]=aTb.data[2];

        minimatrix rblock=minimatrix_blockmatrix_var(this,4,0,3,3);
        minimatrix_memcpy(&rblock,aRb);
#endif // USE_QUATERNIONS


    }



    /// Named constructor with derivatives
    static EssentialMatrix FromRotationAndDirection(const Rot3& aRb, const Unit3& aTb,
            minimatrix* H1,//5, 3
            minimatrix* H2);//5, 2

    /// Named constructor converting a Pose3 with scale to EssentialMatrix (no scale)
    static EssentialMatrix FromPose3(const Pose3& _1P2_,
                                     minimatrix* H);//5, 6

    /// Random, using Rot3::Random and Unit3::Random

    virtual ~EssentialMatrix() {}

    /// @}



    /// @name Manifold
    /// @{
    inline size_t dim() const
    {
        return 5;
    }


    /// Retract delta to manifold
    virtual minimatrix* Retract(const minimatrix* mpose)
    {
        minivector rr(3);
        rr.data[0]=minimatrix_get(mpose,0,0);
        rr.data[1]=minimatrix_get(mpose,1,0);
        rr.data[2]=minimatrix_get(mpose,2,0);

        minivector rt(2);
        rt.data[0]=minimatrix_get(mpose,3,0);
        rt.data[1]=minimatrix_get(mpose,4,0);

        Rot3 rot(rotation());
        Unit3 direct(direction());

        Rot3 rotm;
        minimatrix* rotmr=rot.Retract(&rr);
        minimatrix_memcpy(&rotm,*rotmr);
        delete rotmr;
        Unit3 unittt;
        minimatrix* unitt=direct.Retract(&rt);
        minimatrix_memcpy(&unittt,*unitt);
        delete unitt;

        return new EssentialMatrix(rotm,unittt);
    }

    virtual minimatrix LocalCoordinates(const minimatrix* mpose) const
    {

        EssentialMatrix re;
        minimatrix_memcpy(&re,*mpose);
        Rot3 rot(rotation());
        Rot3 rerot(re.rotation());

        minimatrix v1=rot.LocalCoordinates(&rerot);
        Unit3 redirect(re.direction());
        Unit3 direct(direction());
        minimatrix v2=direct.LocalCoordinates(&redirect);
        minivector result(5);
        result.data[0]=minimatrix_get(v1,0,0);
        result.data[1]=minimatrix_get(v1,1,0);
        result.data[2]=minimatrix_get(v1,2,0);
        result.data[3]=minimatrix_get(v2,0,0);
        result.data[4]=minimatrix_get(v2,1,0);
        return result;
    }
    /// @}

    /// @name Essential matrix methods
    /// @{

    /// Rotation
    inline  Rot3 rotation() const
    {
#ifdef USE_QUATERNIONS
        minivector subq=minimatrix_row(*this,4);
        return Rot3(&subq);
#else
        minimatrix subrot=minimatrix_blockmatrix(*this,4,0,3,3);
        return Rot3(&subrot);
#endif
    }

    /// Direction
    inline Unit3 direction() const
    {
#ifdef USE_QUATERNIONS
        minivector subu=minimatrix_subrow(*this,4,0,3);
        return Unit3(&subu);
#else
        minivector subu=minimatrix_row(*this,4);
        return Unit3(&subu);
#endif // USE_QUATERNIONS
    }

    /// Return 3*3 matrix representation
    inline minimatrix matrix() const
    {
        //return E_;
        //return minimatrix_blockmatrix(*this,0,0,3,3);
        minimatrix result= minimatrix_blockmatrix(*this,0,0,3,3);
        return result;
    }

    /// Return epipole in image_a , as Unit3 to allow for infinity
    inline Unit3 epipole_a() const
    {
#ifdef USE_QUATERNIONS
        minivector subu=minimatrix_subrow(*this,4,0,3);
        return Unit3(&subu);
#else
        minivector subu=minimatrix_row(*this,4);
        return Unit3(&subu);
#endif // USE_QUATERNIONS
    }

    /// Return epipole in image_b, as Unit3 to allow for infinity
    inline Unit3 epipole_b() const
    {
        Unit3 direct(direction());
        return rotation().unrotateUnit(direct);
    }

    /**
     * @brief takes point in world coordinates and transforms it to pose with |t|==1
     * @param p point in world coordinates
     * @param DE optional 3*5 Jacobian wrpt to E
     * @param Dpoint optional 3*3 Jacobian wrpt point
     * @return point in pose coordinates
     */
    minivector transform_to(const minivector& p,
                            minimatrix* DE=NULL,//3, 5
                            minimatrix* Dpoint=NULL) const;//3, 3

    /**
     * Given essential matrix E in camera frame B, convert to body frame C
     * @param cRb rotation from body frame to camera frame
     * @param E essential matrix E in camera frame C
     */
    EssentialMatrix rotate(const Rot3& cRb,minimatrix* HE=NULL, //5, 5
                           minimatrix* HR=NULL) const;// 5, 3

    /**
     * Given essential matrix E in camera frame B, convert to body frame C
     * @param cRb rotation from body frame to camera frame
     * @param E essential matrix E in camera frame C
     */
    friend EssentialMatrix operator*(const Rot3& cRb, const EssentialMatrix& E)
    {
        return E.rotate(cRb);
    }

    /// epipolar error, algebraic
    double error(const minivector& vA,const minivector& vB,
                 minimatrix* H) const;//1, 5

    /// @}



};

};

