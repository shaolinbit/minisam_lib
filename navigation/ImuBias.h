#ifndef IMUBIAS_H
#define IMUBIAS_H


/**
 * @file ImuBias.h
 */

#pragma once
#include <iosfwd>
#include "../mat/Matrix.h"
#include "../mat/MatCal.h"
namespace minisam
{

class ConstantBias:public minimatrix
{

public:

    ConstantBias();
    ~ConstantBias()
    {

    }

    ConstantBias(const  minivector& biasAcc, const  minivector& biasGyro);

    ConstantBias(const ConstantBias& bias);

    ConstantBias(const  minimatrix& v);

    ConstantBias(const minimatrix* v);


    explicit ConstantBias(const  minivector& v);

    ConstantBias(double ax,double ay,double az,double gx,double gy,double gz);

    /** return the accelerometer and gyro biases in a single vector */
    minivector vector() const;

    /** get accelerometer bias */
     minivector accelerometer() const;

    /** get gyroscope bias */
     minivector gyroscope() const;
    /** Correct an accelerometer measurement using this bias model, and optionally compute Jacobians */
    minivector correctAccelerometer(const minivector& measurement,
                                    minimatrix* H1=NULL,
                                    minimatrix* H2=NULL) const;

    /** Correct a gyroscope measurement using this bias model, and optionally compute Jacobians */
    minivector correctGyroscope(const minivector& measurement,
                                minimatrix* H1=NULL,
                                minimatrix* H2=NULL) const;
    /// @}

    /// @name Group
    /// @{

    /** identity for group operation */
    static ConstantBias identity();

    /** inverse */
    ConstantBias operator-() const;

    /** addition of vector on right */
    ConstantBias operator+(const minivector& v) const;

    /** addition */
    ConstantBias operator+(const ConstantBias& b) const;

    /** subtraction */
    ConstantBias operator-(const ConstantBias& b) const;

    /**equal**/
    ConstantBias operator=(const ConstantBias& b) const;

    /// @}

    /// @name Deprecated
    /// @{
    ConstantBias inverse();
    ConstantBias* compose(const ConstantBias& q);
    ConstantBias between(const ConstantBias& q) const;


    virtual minimatrix between(const minimatrix* q) const;
    virtual minimatrix between(const minimatrix* q,minimatrix& H1,minimatrix& H2) const;
    virtual minimatrix* Retract(const minimatrix* mpose);
    virtual minimatrix LocalCoordinates(const minimatrix* mpose) const;


    static minivector Logmap(const ConstantBias& p);
    static ConstantBias Expmap(const minivector& v);
    /// @}


}; // ConstantBias class
ConstantBias ConstantBiasBetween(const ConstantBias& q1,const ConstantBias& q2,
                                 minimatrix* H1,minimatrix* H2);

};
#endif // IMUBIAS_H
