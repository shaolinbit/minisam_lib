/**
 *  @file   GnssBetweenFactor.h
 *  @author Ryan Watson & Jason Gross
 *  @brief  Header file for GnssBetweenFactor
 **/

#pragma once


#include "../gnssNavigation/nonBiasStates.h"
#include "minisam/nonlinear/NonlinearFactor.h"

namespace minisam
{

class  GnssBetweenFactor: public NoiseModelFactor2
{

public:


    GnssBetweenFactor(int state1, int state2, GaussianNoiseModel* model):
        NoiseModelFactor2(model, state1, state2) { }

    virtual NoiseModelFactor* clone() const
    {
        GnssBetweenFactor *newfactor = new GnssBetweenFactor(key1(),key2(),noiseModel());
        return newfactor;
    }


    virtual minivector
    evaluateError(const minimatrix* X1, const minimatrix* X2) const;

    virtual minivector
    evaluateError(const minimatrix* X1, const minimatrix* X2, minimatrix &H1, minimatrix &H2) const;


}; // PseudorangeFactor Factor
}; // namespace
