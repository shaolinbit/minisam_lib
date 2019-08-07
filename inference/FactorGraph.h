#ifndef FACTORPOINTERGRAPH_H
#define FACTORPOINTERGRAPH_H


/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    FactorGraph.h
 * @brief   Factor Graph Base Class
 * @author  Carlos Nieto
 * @author  Christian Potthast
 * @author  Michael Kaess
 * @author  Richard Roberts
 */

#include "../inference/Factor.h"
#include <set>
#include <map>
#include <iostream>
#include <string>


namespace minisam
{
#define VariableSlotsEmpty std::numeric_limits<size_t>::max()
/** Helper */

/**
 * A factor graph is a bipartite graph with factor nodes connected to variable nodes.
 * In this class, however, only factor nodes are kept around.
 * \nosubgrouping
 */

template<typename TPFactor>
class FactorGraph
{

public:
    std::vector<TPFactor*> factors_;
    typedef typename std::vector<TPFactor*>::iterator iterator;
    typedef typename std::vector<TPFactor*>::const_iterator const_iterator;
public:
    /** Default constructor */
    FactorGraph() {}

    FactorGraph(const FactorGraph&) = default;

    /** Constructor from iterator over factors (shared_ptr or plain objects) */
    FactorGraph(iterator firstFactor, iterator lastFactor)
    {
        push_back(firstFactor, lastFactor);
    }

    /** Construct from container of factors (shared_ptr or plain objects) */
    template<class CONTAINER>
    explicit FactorGraph(const CONTAINER& factors)
    {
        for(const_iterator pfactor = factors.begin(); pfactor != factors.end(); ++pfactor)
        {
            factors_.push_back(pfactor);
        }
    }

    FactorGraph(const std::vector<TPFactor*> factors)
    {
        for(const_iterator pfactor = factors.begin(); pfactor != factors.end(); ++pfactor)
        {
            factors_.push_back(*pfactor);
        }
    }

    FactorGraph(FactorGraph& FG):FactorGraph(FG.factors_) {}

    /// @}
    /// @name Advanced Constructors
    /// @{

public:
    /// @name Adding Factors
    /// @{

    /**
     * Reserve space for the specified number of factors if you know in
     * advance how many there will be (works like FastVector::reserve).
     */
    void reserve(int size)
    {
        factors_.reserve(size);
    }

    /** Add a factor directly using a shared_ptr */

    void push_back(const TPFactor* factor)
    {
        factors_.push_back(factor);
    }

    void push_back(TPFactor* factor)
    {
        factors_.push_back(factor);
    }

    /** push back many factors with an iterator over shared_ptr (factors are not copied) */
    void  push_back(iterator firstFactor, iterator lastFactor)
    {
        factors_.insert(end(), firstFactor, lastFactor);
    }


public:
    /// @name Standard Interface
    /// @{

    /** return the number of factors (including any null factors set by remove() ). */
    int size() const
    {
        return factors_.size();
    }

    /** Check if the graph is empty (null factors set by remove() will cause this to return false). */
    bool empty() const
    {
        return factors_.empty();
    }

    /** Get a specific factor by index (this checks array bounds and may throw an exception, as
     *  opposed to operator[] which does not).
     */
    const TPFactor* at(int i) const
    {
        return factors_.at(i);
    }

    /** Get a specific factor by index (this checks array bounds and may throw an exception, as
     *  opposed to operator[] which does not).
     */
    TPFactor* at(int i)
    {
        return factors_.at(i);
    }
    TPFactor* nconstat(int i) const
    {
        return factors_.at(i);
    }

    TPFactor* nconstat(int i)
    {
        return factors_.at(i);
    }

    /** Get a specific factor by index (this does not check array bounds, as opposed to at() which
     *  does).
     */
    const TPFactor* operator[](int i) const
    {
        return at(i);
    }

    /** Get a specific factor by index (this does not check array bounds, as opposed to at() which
     *  does).
     */
    TPFactor* operator[](int i)
    {
        return at(i);
    }

    /** Iterator to beginning of factors. */
    const_iterator begin() const
    {
        return factors_.begin();
    }

    /** Iterator to end of factors. */
    const_iterator end()   const
    {
        return factors_.end();
    }

    /** Get the first factor */
    TPFactor* front() const
    {
        return factors_.front();
    }

    /** Get the last factor */
    TPFactor* back() const
    {
        return factors_.back();
    }

    /// @}
    /// @name Modifying Factor Graphs (imperative, discouraged)
    /// @{

    /** non-const STL-style begin() */
    iterator begin()
    {
        return factors_.begin();
    }

    /** non-const STL-style end() */
    iterator end()
    {
        return factors_.end();
    }

    /** Directly resize the number of factors in the graph. If the new size is less than the
     * original, factors at the end will be removed.  If the new size is larger than the original,
     * null factors will be appended.
     */
    void resize(int size)
    {
        factors_.resize(size);
    }
    void clear()
    {
        factors_.clear();
    }

    /** delete factor without re-arranging indexes by inserting a NULL pointer */
    void remove(int i)
    {
        factors_[i]->keys_.clear();
    }


    /** replace a factor by index */
    void replace(int index, TPFactor* factor)
    {
        at(index) = factor;
    }

    /** Erase factor and rearrange other factors to take up the empty space */
    iterator erase(iterator item)
    {
        return factors_.erase(item);
    }

    iterator insert(iterator ipos,TPFactor* factor)
    {
        return factors_.insert(ipos,factor);
    }

    /** Erase factors and rearrange other factors to take up the empty space */
    iterator erase(iterator first, iterator last)
    {
        return factors_.erase(first, last);
    }

    /// @}
    /// @name Advanced Interface
    /// @{

    /** return the number of non-null factors */
    int nrFactors() const
    {
        int size_ = 0;
        for(const TPFactor* factor: factors_)
            size_++;
        return size_;
    }

    /** Potentially slow function to return all keys involved, sorted, as a set */
    std::set<int> keys() const
    {
        std::set<int> keys;
        for(const TPFactor* factor: this->factors_)
        {
            keys.insert(factor->begin(), factor->end());
        }
        return keys;
    }


    FactorGraph& operator=( FactorGraph& other)
    {
        factors_.clear();
        for(auto& tf:other)
        {
            this->push_back(tf);
        }
    }

    FactorGraph& operator=( const FactorGraph& other)
    {
        factors_.clear();
        for(auto& tf:other)
        {
            this->push_back(tf);
        }
    }
};
};

#endif
