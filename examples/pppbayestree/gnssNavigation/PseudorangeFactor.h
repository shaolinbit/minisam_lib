/**
 *  @file   PseudorangeFactor.h
 *  @author Ryan Watson & Jason Gross
 *  @brief  Header file for Pseudorange factor
 **/

#pragma once

#include "../gnssNavigation/GnssTools.h"
#include "minisam/nonlinear/NonlinearFactor.h"
#include "../gnssNavigation/nonBiasStates.h"

namespace minisam
{

class  PseudorangeFactor : public NoiseModelFactor1
{
private:
    typedef NoiseModelFactor1 Base;
    minivector* nomXYZ_;
    minivector* satXYZ_;
    minivector* h_;
    //nonBiasStates h_;
    double measured_;

public:

    typedef PseudorangeFactor This;

    PseudorangeFactor() : measured_(0),nomXYZ_(new minivector(3)),
    satXYZ_(new minivector(3)),h_(new nonBiasStates(1,1,1,1,1))
    {
       
    }

    virtual ~PseudorangeFactor()
    {

     delete nomXYZ_;
     delete satXYZ_;
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


    PseudorangeFactor(int key, const double deltaObs, minivector* satXYZ,
                      minivector* nomXYZ, GaussianNoiseModel* model) :
        Base(model, key), measured_(deltaObs),h_(new nonBiasStates(1,1,1,1,1)),
        satXYZ_(satXYZ),nomXYZ_(nomXYZ)
    {
       
    }


    virtual NoiseModelFactor* clone() const
    {
        PseudorangeFactor* nprf=new PseudorangeFactor(key(),measured_,satXYZ_,nomXYZ_,noiseModel());
        return nprf;

    }


    virtual minivector  evaluateError(const minimatrix* q) const;


    virtual minivector  evaluateError(const minimatrix* q,
                              minimatrix& H) const;


}; // PseudorangeFactor Factor
}; // namespace
