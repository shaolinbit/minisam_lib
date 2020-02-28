#ifndef BETWEENFACTOR_H
#define BETWEENFACTOR_H

/**
 *  @file  BetweenFactor.h
 **/
#pragma once
#include "../nonlinear/NonlinearFactor.h"
namespace minisam
{
/**
 * A class for a measurement predicted by "between(config[key1],config[key2])"
 * @addtogroup SLAM
 */
class BetweenFactor: public NoiseModelFactor2
{

public:
    minimatrix* measured_;

public:

    /** default constructor - only use for serialization */
    BetweenFactor() {}

    /** Constructor */
    BetweenFactor(int key1, int key2, minimatrix* measured,
                  GaussianNoiseModel* model) :
        NoiseModelFactor2(model, key1, key2), measured_(measured)
    {
    }

    ~BetweenFactor()
    {
        delete measured_;
    }


    virtual minivector evaluateError(const minimatrix* p1, const minimatrix* p2) const
    {
        minimatrix hx =p1->between(p2);
        minimatrix resultm=measured_->LocalCoordinates(&hx);
        return minivector(resultm);
    }
    virtual  minivector evaluateError(const minimatrix* p1,const minimatrix* p2,
                                      minimatrix& H1,minimatrix& H2) const
    {
        minimatrix hx =p1->between(p2,H1,H2);
        minimatrix resultm=measured_->LocalCoordinates(&hx,&H1,&H2);
        return minivector(resultm);
    }

    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x,
                                       std::vector<minimatrix> &H) const
    {
        std::map<int,minimatrix*>::const_iterator itb1=x.find(key1());
        std::map<int,minimatrix*>::const_iterator itb2=x.find(key2());

        return evaluateError(itb1->second,itb2->second,
                             *(H.begin()),*(H.begin()+1));
    }

    virtual minivector unwhitenedError(const std::map<int, minimatrix*>& x) const
    {
        std::map<int,minimatrix*>::const_iterator itb1=x.find(key1());
        std::map<int,minimatrix*>::const_iterator itb2=x.find(key2());

        return evaluateError(itb1->second,itb2->second);
    }


    /** return the measured */
    const minimatrix& measured() const
    {
        return *measured_;
    }

    /** number of variables attached to this factor */
    int size() const
    {
        return 2;
    }

}; // \class BetweenFactor
};

#endif // BETWEENFACTOR_H
