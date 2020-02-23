/**
 *  @file   PhaseFactor.h
 *  @author Ryan Watson & Jason Gross
 *  @brief  Header file for Carrier-Phase Factor
 **/

#pragma once

#include "../gnssNavigation/GnssTools.h"
#include "minisam/nonlinear/NonlinearFactor.h"
#include "../gnssNavigation/nonBiasStates.h"

namespace minisam
{

class  PhaseFactor : public NoiseModelFactor2
{
private:
    typedef NoiseModelFactor2 Base;
    minivector* satXYZ_;
    minivector* nomXYZ_;
    double measured_;
    //nonBiasStates h_;
    minivector* h_;

public:

    typedef PhaseFactor This;

    PhaseFactor() : measured_(0),satXYZ_(new minivector(3)),
    nomXYZ_(new minivector(3)),h_(new nonBiasStates(1,1,1,1,1))
    {
       
    }

    virtual ~PhaseFactor()
    {
        delete satXYZ_;
        delete nomXYZ_;
        delete h_;
        delete noiseModel_;
    }

    minivector* satXYZ()
    {
        return satXYZ_;
    }

    minivector* nomXYZ()
    {
        return nomXYZ_;
    }
    double measured()
    {
        return measured_;
    }

    PhaseFactor(int deltaStates, int bias, const double measurement,
                minivector* satXYZ, minivector* nomXYZ,GaussianNoiseModel* model) :
        Base(model, deltaStates, bias), measured_(measurement),satXYZ_(satXYZ),nomXYZ_(nomXYZ),
        h_(new nonBiasStates(1,1,1,1,1))
    {
       
    }


    virtual NoiseModelFactor* clone() const
    {
        PhaseFactor* npf=new PhaseFactor(key1(),key2(),measured_,satXYZ_,nomXYZ_,noiseModel());
        return npf;
    }

    virtual minivector evaluateError(const minimatrix* q, const minimatrix* g) const;

    virtual minivector evaluateError(const minimatrix* q, const minimatrix* g,
                             minimatrix& H1,minimatrix& H2) const;


}; // PhaseFactor Factor
}; // namespace
