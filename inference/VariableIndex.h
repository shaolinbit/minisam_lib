#ifndef VARIABLEINDEX_H
#define VARIABLEINDEX_H


/**
 * @file    VariableIndex.h
 * @author
 * @date
 */

//#pragma once

#include "../inference/Factor.h"
//#include "../inference/FactorGraph.h"
//#include "../linear/GaussianFactorGraph.h"
#include "../linear/GaussianFactorGraph.h"
#include "../nonlinear/NonlinearFactor.h"
//#include "../nonlinear/NonlinearFactorGraph.h"
#include "../nonlinear/NonlinearFactorGraph.h"
//#include "../nonlinear/NoiseModelFactorGraph.h"
#include <map>
#include <cassert>


/**
 * The VariableIndex class computes and stores the block column structure of a
 * factor graph.  The factor graph stores a collection of factors, each of
 * which involves a set of variables.  In contrast, the VariableIndex is built
 * from a factor graph prior to elimination, and stores the list of factors
 * that involve each variable.  This information is stored as a deque of
 * lists of factor indices.
 * \nosubgrouping
 */
//class GaussianFactorGraph;
class GaussianFactorGraph;
class NonlinearFactorGraph;

class  VariableIndex
{
public:
    std::map<int,std::vector<int>> index_;
public:
    int nFactors_; // Number of factors in the original factor graph.
    int nEntries_; // Sum of involved variable counts of each factor.
public:

    /// @name Standard Constructors
    /// @{

    /** Default constructor, creates an empty VariableIndex */
    VariableIndex() : nFactors_(0), nEntries_(0) {}

    ~VariableIndex() {}

    VariableIndex(const VariableIndex& other)
    {
        *this=other;
    }

    /**
     * Create a VariableIndex that computes and stores the block column structure
     * of a factor graph.
     */
    //template<class TPFactorGraph>
   // VariableIndex(const TPFactorGraph& factorGraph) : nFactors_(0), nEntries_(0)
   /// {
  //      augment(factorGraph);
  //  }


    /*
    VariableIndex(const GaussianFactorGraph& factorGraph) : nFactors_(0), nEntries_(0)
    {
        augment(factorGraph);
    }*/

    /*
    VariableIndex(const FactorGraph<NonlinearFactor>& factorGraph) : nFactors_(0), nEntries_(0)
    {
        augment(factorGraph);
    }*/

    VariableIndex(const NonlinearFactorGraph& factorGraph) : nFactors_(0), nEntries_(0)
    {
        augment(factorGraph);
    }

    VariableIndex(const FactorPointerGraph<NonlinearFactor>& factorGraph) : nFactors_(0), nEntries_(0)
    {
        augment(factorGraph);
    }
    VariableIndex(const GaussianFactorGraph& factorGraph) : nFactors_(0), nEntries_(0)
    {
        augment(factorGraph);
    }

    /// @}
    /// @name Standard Interface
    /// @{

    /**
     * The number of variable entries.  This is one greater than the variable
     * with the highest index.
     */
    int size() const
    {
        return index_.size();
    }

    /** The number of factors in the original factor graph */
    int nFactors() const
    {
        return nFactors_;
    }

    /** The number of nonzero blocks, i.e. the number of variable-factor entries */
    int nEntries() const
    {
        return nEntries_;
    }

    /** Access a list of factors by variable */

    const std::vector<int>& operator[](int variable) const
    {
        std::map<int,std::vector<int>>::const_iterator item = index_.find(variable);
        if(item == index_.end())
            throw std::invalid_argument("Requested non-existent variable from VariableIndex");
        else
            return item->second;
        //return item->second;
    }

    VariableIndex& operator=(const VariableIndex& other)
    {
        index_=other.index_;
        nFactors_=other.nFactors_;
        nEntries_=other.nEntries_;

        return *this;
    }

    /// @}
    /// @name Testable
    /// @{

    /** Test for equality (for unit tests and debug assertions). */
// bool equals(const VariableIndex& other, double tol=0.0) const;

    /** Print the variable index (for unit tests and debugging). */
    void print(const std::string& str = "VariableIndex: ") const;

    /**
     * Output dual hypergraph to Metis file format for use with hmetis
     * In the dual graph, variables are hyperedges, factors are nodes.
     */
    void outputMetisFormat(std::ostream& os) const;


// const std::vector<int>& operator[](int j) const
// {
//     return index_[j];
//  }


    /// @}
    /// @name Advanced Interface
    /// @{

    /**
     * Augment the variable index with new factors.  This can be used when
     * solving problems incrementally.
     */
   // template<class TPFactorGraph>
   // void augment(const TPFactorGraph& factors, const std::vector<int>& newFactorIndices);

    void augmentnfpindex(NonlinearFactorGraph& factors,std::vector<int>& newFactorIndices);

   // template<class TPFactorGraph>
   // void augment( TPFactorGraph& factors, std::vector<int>& newFactorIndices);

   // template<class TPFactorGraph>
   // void augment(const TPFactorGraph& factors);

   // template<class TPFactorGraph>
  //  void augment(TPFactorGraph& factors);

   // void augment(const GaussianFactorGraph& factors);
    void augment(const GaussianFactorGraph& factors);

   // void augment(const FactorGraph<NonlinearFactor>& factors);
  //  void augment(const FactorGraph<NoiseModelFactor>& factors);

    //void augment(const FactorGraph<RealGaussianFactor>& factors);

    void augment(const FactorPointerGraph<NonlinearFactor>& factors);

    void augment(const GaussianFactorGraph* factors);

//void augment(const NonlinearFactorGraph& factors);

    //void augment(const NoiseModelFactorGraph& factors);

    void augment(const NonlinearFactorGraph& factors);

    void augmentnfp(NonlinearFactorGraph& factors);

    /**
     * Remove entries corresponding to the specified factors. NOTE: We intentionally do not decrement
     * nFactors_ because the factor indices need to remain consistent.  Removing factors from a factor
     * graph does not shift the indices of other factors.  Also, we keep nFactors_ one greater than
     * the highest-numbered factor referenced in a VariableIndex.
     *
     * @param indices The indices of the factors to remove, which must match \c factors
     * @param factors The factors being removed, which must symbolically correspond exactly to the
     *        factors with the specified \c indices that were added.
     */
   // template<class TPFactorGraph>
   // void remove(std::vector<int>::iterator firstFactor, std::vector<int>::iterator lastFactor, const TPFactorGraph& factors);

   // template<class TPFactorGraph>
   // void remove(std::set<int>::iterator firstFactor, std::set<int>::iterator lastFactor, const TPFactorGraph& factors);

   // template<class TPFactorGraph>
   // void remove(std::set<int>::iterator firstFactor, std::set<int>::iterator lastFactor, TPFactorGraph& factors);


   // template<class TPFactorGraph>
   // void remove(std::vector<int>::const_iterator firstFactor, std::vector<int>::const_iterator lastFactor, const TPFactorGraph& factors);

    //template<class TPFactorGraph>
   // void remove(std::vector<int>::iterator firstFactor, std::vector<int>::iterator lastFactor, TPFactorGraph& factors);

    void remove(std::set<int>::iterator firstFactor, std::set<int>::iterator lastFactor, FactorPointerGraph<Factor>& factors);

    void removenf(std::vector<int>::iterator firstFactor, std::vector<int>::iterator lastFactor, const NonlinearFactorGraph& factors);


    //template<class TPFactorGraph>
    //void removenl(std::vector<int>::iterator firstFactor, std::vector<int>::iterator lastFactor, NonlinearFactorGraph& factors);


    /** Remove unused empty variables (in debug mode verifies they are empty). */
    //template<typename ITERATOR>
    void removeUnusedVariables(std::set<int>::iterator firstKey, std::set<int>::iterator lastKey);

    void newremoveUnusedVariables(const std::set<int>& unusedkey);


    /** Iterator to the first variable entry */
    std::map<int,std::vector<int>>::const_iterator begin() const
    {
        return index_.begin();
    }

    /** Iterator to the first variable entry */
    std::map<int,std::vector<int>>::const_iterator end() const
    {
        return index_.end();
    }

    /** Find the iterator for the requested variable entry */
    std::map<int,std::vector<int>>::const_iterator find(int key) const
    {
        return index_.find(key);
    }

protected:
    std::vector<int>::iterator factorsBegin(int variable)
    {
        return internalAtnc(variable).begin();
    }
    std::vector<int>::iterator factorsEnd(int variable)
    {
        return internalAtnc(variable).end();
    }

    std::vector<int>::const_iterator factorsBegin(int variable) const
    {
        return internalAt(variable).begin();
    }
    std::vector<int>::const_iterator factorsEnd(int variable) const
    {
        return internalAt(variable).end();
    }

    /// Internal version of 'at' that asserts existence
    const std::vector<int>& internalAt(int variable) const
    {
        const std::map<int,std::vector<int>>::const_iterator item = index_.find(variable);
        assert(item != index_.end());
        return item->second;
    }


    /// Internal version of 'at' that asserts existence
    std::vector<int>& internalAtnc(int variable)
    {
        const std::map<int,std::vector<int>>::iterator item = index_.find(variable);
        assert(item != index_.end());
        return item->second;
    }

    /// @}
};
#endif // VARIABLEINDEX_H
