#ifndef SYMBOLICCONDITONAL_H
#define SYMBOLICCONDITONAL_H

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    SymbolicConditional.h
 * @author  Richard Roberts
 * @date    Oct 17, 2010
 */


#include "../inference/Factor.h"

namespace minisam
{

/**
 * SymbolicConditional is a conditional with keys but no probability
 * data, produced by symbolic elimination of SymbolicFactor.
 * \nosubgrouping
 */
class  SymbolicConditional :public Factor
{
public:
    Factor* nrFrontals_;
    Factor* nrParents_;
public:

    /// @name Standard Constructors
    /// @{

    /** Empty Constructor to make serialization possible */
    SymbolicConditional();// {}

    /** No parents */
    SymbolicConditional(const Factor& fFactor) : Factor(fFactor), nrFrontals_(new Factor(fFactor)),nrParents_(new Factor()) {}

    SymbolicConditional(const Factor& fFactor,int sizeFrontal);

    SymbolicConditional(const SymbolicConditional& other);

    ~SymbolicConditional();// {}
    int getsizenrFrontals() const
    {
        return nrFrontals_->size();
    }

    int getsizenrParents() const
    {
        return nrParents_->size();
    }

    std::vector<int>::const_iterator cbeginFrontals()
    const;

    std::vector<int>::const_iterator cendFrontals()
    const;

    std::vector<int>::const_iterator cbeginParents()
    const;

    std::vector<int>::const_iterator cendParents()
    const;

    /// @}
    /// @name Advanced Interface
    /// @{

    Factor* nrFrontals() const;
    Factor* nrParents() const;


    /// @}
};
};
#endif // SYMBOLICCONDITONAL_H
