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
    int Frontalsize_;
    int Parentsize_;
public:

    /// @name Standard Constructors
    /// @{

    /** Empty Constructor to make serialization possible */
    SymbolicConditional();// {}

    /** No parents */
    SymbolicConditional(const Factor& fFactor);// : Factor(fFactor), nrFrontals_(new Factor(fFactor)),nrParents_(new Factor()) {}

    SymbolicConditional(const Factor& fFactor,int sizeFrontal);

    SymbolicConditional(const SymbolicConditional& other);

    ~SymbolicConditional();

    /// @}
};
};
#endif // SYMBOLICCONDITONAL_H
