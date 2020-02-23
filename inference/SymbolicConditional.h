#ifndef SYMBOLICCONDITONAL_H
#define SYMBOLICCONDITONAL_H

/**
 * @file    SymbolicConditional.h
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
    SymbolicConditional();

    /** No parents */
    SymbolicConditional(const Factor& fFactor);

    SymbolicConditional(const Factor& fFactor,int sizeFrontal);

    SymbolicConditional(const SymbolicConditional& other);

    ~SymbolicConditional();

    /// @}
};
};
#endif // SYMBOLICCONDITONAL_H
