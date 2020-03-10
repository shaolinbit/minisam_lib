#ifndef PRIORFACTOR_H
#define PRIORFACTOR_H

/**
 *  @file  PriorFactor.h
 **/
#pragma once

#include "../nonlinear/NonlinearFactor.h"
#include "../linear/NoiseModel.h"
namespace minisam
{

/**
 * A class for a soft prior
 * @addtogroup SLAM
 */

class PriorFactor: public NoiseModelFactor1
{

public:
    minimatrix* prior_; /** The measurement */
public:
    /** default constructor - only use for serialization */
    PriorFactor() {}

    virtual ~PriorFactor()
    {
       if(prior_!=NULL)
        {
        delete prior_;
        prior_=NULL;
        }
    }

    /** Constructor */
    PriorFactor(int key, minimatrix* prior, GaussianNoiseModel* model) :
        NoiseModelFactor1(model, key), prior_(prior)
    {
    }

    /** Convenience constructor that takes a full covariance argument */
    PriorFactor(int key, minimatrix* prior, const minimatrix& covariance) :
        NoiseModelFactor1(GaussianNoiseModel_Covariance(covariance), key), prior_(prior)
    {}


    /** implement functions needed for Testable */

    virtual minivector unwhitenedError(const std::map<int,minimatrix*>& x)const
    {
        std::map<int,minimatrix*>::const_iterator itb=x.find(key());
        return evaluateError(itb->second);
    }
    virtual minivector unwhitenedError(const std::map<int,minimatrix*>& x,std::vector<minimatrix>& H) const
    {
        std::map<int,minimatrix*>::const_iterator itb=x.find(key());
        return evaluateError(itb->second,*(H.begin()));
    }

    virtual minivector evaluateError(const minimatrix* x) const
    {
        minimatrix v=x->LocalCoordinates(prior_);
        minimatrix_scale(&v,-1.0);
        minivector result(v);
        return result;
    }
    virtual minivector evaluateError(const minimatrix* x, minimatrix& H) const
    {
        minimatrix_resize(&H,x->dimension,x->dimension);
        minimatrix_set_identity(&H);

        minimatrix v=x->LocalCoordinates(prior_);
        minimatrix_scale(&v,-1.0);
        minivector result(v);

        if(DEBUGSTATE)
        {
        cout<<"key()"<<endl;
        cout<<key()<<endl;
        cout<<"x"<<endl;
        minimatrix_print(x);
        cout<<"prior_"<<endl;
        minimatrix_print(prior_);
        }

        return result;
    }




    const minimatrix* prior() const
    {
        return prior_;
    }


    virtual NoiseModelFactor* clone()const
    {
        PriorFactor* newfactor=new PriorFactor(key(),prior_,noiseModel_);
        return newfactor;
    }
};
};

#endif
