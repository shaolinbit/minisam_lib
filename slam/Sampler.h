#ifndef SAMPLER_H
#define SAMPLER_H


/**
 * @brief sampling that can be parameterized using a NoiseModel to generate samples from
 * @file Sampler.h
 * the given distribution
 */

#pragma once

#include "../linear/NoiseModel.h"


namespace minisam {

/**
 * Sampling structure that keeps internal random number generators for
 * diagonal distributions specified by NoiseModel
 *
 * This is primarily to allow for variable seeds, and does roughly the same
 * thing as sample() in NoiseModel.
 */
class  Sampler
{
protected:
    /** noiseModel created at generation */
    GaussianNoiseModel* model_;

    /** generator */
    int generatoridium_;

public:

    /**
     * Create a sampler for the distribution specified by a diagonal NoiseModel
     * with a manually specified seed
     *
     * NOTE: do not use zero as a seed, it will break the generator
     */
    Sampler(GaussianNoiseModel* model, int seed = 42u);

    /**
     * Create a sampler for a distribution specified by a vector of sigmas directly
     *
     * NOTE: do not use zero as a seed, it will break the generator
     */
    Sampler(const minivector& sigmas, int seed = 42u);

    /**
     * Create a sampler without a given noisemodel - pass in to sample
     *
     * NOTE: do not use zero as a seed, it will break the generator
     */
    Sampler(int seed = 42u);

    /** access functions */
    int dim() const
    {
        //assert(model_.get());
        return model_->dim();
    }
    minivector sigmas() const
    {
        // assert(model_.get());
        return model_->sigmas();
    }
    GaussianNoiseModel* model() const
    {
        return model_;
    }

    /**
     * sample from distribution
     * NOTE: not const due to need to update the underlying generator
     */
    minivector  sample();

    /**
     * Sample from noisemodel passed in as an argument,
     * can be used without having initialized a model for the system.
     *
     * NOTE: not const due to need to update the underlying generator
     */
    minivector  sampleNewModel(GaussianNoiseModel* model);

protected:

    /** given sigmas for a diagonal model, returns a sample */
    minivector  sampleDiagonal(const minivector & sigmas);

};
double ran1(int *idum);
double UnitWhiteRand(int *idum);
double NormalRand(int *idum,double mean,double cov);


}; // \namespace minisam
#endif // SAMPLER_H
