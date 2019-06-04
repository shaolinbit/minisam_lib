#ifndef FACTORPOINTERGRAPH_H
#define FACTORPOINTERGRAPH_H

/**
 * @file    FactorPointerGraph.h
 * @brief   Factor Graph Base Class
 * @author
 */

// \callgraph

//#pragma once

#include "../inference/Factor.h"
#include <set>
#include <map>
#include <iostream>
#include <string>

#define VariableSlotsEmpty std::numeric_limits<size_t>::max()
/** Helper */

/**
 * A factor graph is a bipartite graph with factor nodes connected to variable nodes.
 * In this class, however, only factor nodes are kept around.
 * \nosubgrouping
 */
template<typename TPFactor>
class FactorPointerGraph
{

public:
    //  typedef FACTOR FactorType;  ///< factor type
    //typedef boost::shared_ptr<FACTOR> sharedFactor;  ///< Shared pointer to a factor
    //typedef sharedFactor value_type;
    // typedef typename FastVector<sharedFactor>::iterator iterator;
    // typedef typename FastVector<sharedFactor>::const_iterator const_iterator;
    std::vector<TPFactor*> factors_;
    typedef typename std::vector<TPFactor*>::iterator iterator;
    typedef typename std::vector<TPFactor*>::const_iterator const_iterator;
public:
    /** Default constructor */
    FactorPointerGraph() {}

    FactorPointerGraph(const FactorPointerGraph&) = default;

    /** Constructor from iterator over factors (shared_ptr or plain objects) */
    //template<class TPFactor>
    FactorPointerGraph(iterator firstFactor, iterator lastFactor)
    {
        push_back(firstFactor, lastFactor);
    }

    /** Construct from container of factors (shared_ptr or plain objects) */
    template<class CONTAINER>
    explicit FactorPointerGraph(const CONTAINER& factors)
    {
        for(const_iterator pfactor = factors.begin(); pfactor != factors.end(); ++pfactor)
        {
            factors_.push_back(pfactor);
        }
    }

    FactorPointerGraph(const std::vector<TPFactor*> factors)
    {
        for(const_iterator pfactor = factors.begin(); pfactor != factors.end(); ++pfactor)
        {
            factors_.push_back(*pfactor);
        }
    }

    FactorPointerGraph(FactorPointerGraph& FG):FactorPointerGraph(FG.factors_) {}

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

    // TODO: are these needed?

    /** Add a factor directly using a shared_ptr */
    //template<class DER GaussianFactorGraph clone() const;IVEDFACTOR>
    // typename std::enable_if<std::is_base_of<FactorType, DERIVEDFACTOR>::value>::type
    // push_back(Factor* factor) {
    //  factors_.push_back(factor); }

    /** Add a factor directly using a shared_ptr */

    void push_back(const TPFactor* factor)
    {
        factors_.push_back(factor);
    }

    void push_back(TPFactor* factor)
    {
        factors_.push_back(factor);
    }


    /** Emplace a factor */
    // template<class DERIVEDFACTOR, class... Args>
    //  typename std::enable_if<std::is_base_of<FactorType, DERIVEDFACTOR>::value>::type
    //  void emplace_shared(Args&&... args) {
    //      factors_.push_back(Factor*(std::forward<Args>(args)...));
//   }

    /** push back many factors with an iterator over shared_ptr (factors are not copied) */
    //template<typename ITERATOR>
    //typename std::enable_if<std::is_base_of<FactorType, typename ITERATOR::value_type::element_type>::value>::type
    void  push_back(iterator firstFactor, iterator lastFactor)
    {
        factors_.insert(end(), firstFactor, lastFactor);
    }

    /** push back many factors as shared_ptr's in a container (factors are not copied) */
    // template<typename CONTAINER>
    // typename std::enable_if<std::is_base_of<FactorType, typename CONTAINER::value_type::element_type>::value>::type
    //  void push_back(const std::vector<Factor>& container) {
    //    push_back(container.begin(), container.end());
    // }

    /** push back a BayesTree as a collection of factors.  NOTE: This should be hidden in derived
     *  classes in favor of a type-specialized version that calls this templated function. */
    // template<class CLIQUE>
    // typename std::enable_if<std::is_base_of<This, typename CLIQUE::FactorGraphType>::value>::type


    //wait for the definition of bayestree
    // void  push_back_factor(const BayesTree<CLIQUE>& bayesTree)
    //  {
    //    bayesTree.addFactorsToGraph(*this);
    //  }


    /** push back many factors with an iterator over plain factors (factors are copied) */
    // template<typename ITERATOR>
    // typename std::enable_if<std::is_base_of<FactorType, typename ITERATOR::value_type>::value>::type
    // void   push_back_factor(std::vector<Factor>::iterator firstFactor, std::vector<Factor>::iterator  lastFactor) {
    //      for ((std::vector<Factor>::iterator f = firstFactor; f != lastFactor; ++f)
    //        push_back_factor(*f);
    //  }
//
    /** push back many factors as non-pointer objects in a container (factors are copied) */
    //template<typename CONTAINER>
    //typename std::enable_if<std::is_base_of<FactorType, typename CONTAINER::value_type>::value>::type
    //void  push_back_factor(const std::vector<Factor>& container) {
    //    push_back_factor(container.begin(), container.end());
    //}

    /** Add a factor directly using a shared_ptr */
    //  template<class DERIVEDFACTOR>
    /*
      typename std::enable_if<std::is_base_of<FactorType, DERIVEDFACTOR>::value,
        boost::assign::list_inserter<RefCallPushBack<This> > >::type
        operator+=(boost::shared_ptr<DERIVEDFACTOR> factor) {
          return boost::assign::make_list_inserter(RefCallPushBack<This>(*this))(factor);
      }


      boost::assign::list_inserter<CRefCallPushBack<This> >
        operator+=(const sharedFactor& factor) {
          return boost::assign::make_list_inserter(CRefCallPushBack<This>(*this))(factor);
      }


      template<class FACTOR_OR_CONTAINER>
      boost::assign::list_inserter<CRefCallPushBack<This> >
        operator+=(const FACTOR_OR_CONTAINER& factorOrContainer) {
          return boost::assign::make_list_inserter(CRefCallPushBack<This>(*this))(factorOrContainer);
      }
      */

    /** Add a factor directly using a shared_ptr */
    // template<class DERIVEDFACTOR>
    // typename std::enable_if<std::is_base_of<FactorType, DERIVEDFACTOR>::value>::type
    void  add(std::vector<TPFactor*>& factors)
    {
        //factors_.push_back(factor);

        for(const_iterator pfactor = factors.begin(); pfactor != factors.end(); ++pfactor)
        {
            //cout << *iter << endl;
            factors_.push_back(*pfactor);
        }
    }

    /** Add a factor directly using a shared_ptr */
    void add(const std::vector<TPFactor*>& factor)
    {
        //for_each(*factor.begin(),*factor.end(),factors_.push_back  )
        // push_back_factor(factor);
        for(const_iterator pfactor = factor.begin(); pfactor != factor.end(); ++pfactor)
        {
            //cout << *iter << endl;
            factors_.push_back(*pfactor);
        }
    }

    /** Add a factor or container of factors, including STL collections, BayesTrees, etc. */
    //template<class FACTOR_OR_CONTAINER>
    //void add(const FACTOR_OR_CONTAINER& factorOrContainer) {
    //  push_back_factor(factorOrContainer);
    //  }

    /// @}
    /// @name Testable
    /// @{

    /** print out graph */
    // void print(const std::string& s = "FactorGraph",
    //   const KeyFormatter& formatter = DefaultKeyFormatter) const;

    /** Check equality */
    //  bool equals(const This& fg, double tol = 1e-9) const;
    /// @}

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
        // factors_[i].reset();
        // factors_.erase(factors_.begin()+i);
        //factors_[i]=NULL;
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

    /*
    void insert (iterator position, TPFactor* val)
    {
        factors_.insert(position,1,val);
    }*/
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

    /** Potentially slow function to return all keys involved, sorted, as a vector
    std::vector<int> keyVector() const
    {
        std::vector<int> keys;
        keys.reserve(2 * size());  // guess at size
        for (auto factor: factors_)
            keys.insert(keys.end(), factor->begin(), factor->end());
        std::sort(keys.begin(), keys.end());
        auto last = std::unique(keys.begin(), keys.end());
        keys.erase(last, keys.end());
        return keys;
    }*/

    FactorPointerGraph& operator=( FactorPointerGraph& other)
    {
        factors_.clear();
        for(auto& tf:other)
        {
            this->push_back(tf);
        }
    }

    FactorPointerGraph& operator=( const FactorPointerGraph& other)
    {
        factors_.clear();
        for(auto& tf:other)
        {
            this->push_back(tf);
        }
    }
};


#endif
