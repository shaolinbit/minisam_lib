#ifndef ELIMINATIONTREEPOINTER_H_INCLUDED
#define ELIMINATIONTREEPOINTER_H_INCLUDED

/**
* @file    EliminationTree.h
* @author
* @date
*/
//#pragma once

#include <utility>
#include "../linear/GaussianFactorGraph.h"
#include "../linear/RealGaussianFactor.h"
#include "../linear/GaussianBayesNetPointer.h"

class VariableIndex;
class Ordering;
class GaussianFactorGraph;

/**
* An elimination tree is a data structure used intermediately during
* elimination.  In future versions it will be used to save work between
* multiple eliminations.
*
* When a variable is eliminated, a new factor is created by combining that
* variable's neighboring factors.  The new combined factor involves the combined
* factors' involved variables.  When the lowest-ordered one of those variables
* is eliminated, it consumes that combined factor.  In the elimination tree,
* that lowest-ordered variable is the parent of the variable that was eliminated to
* produce the combined factor.  This yields a tree in general, and not a chain
* because of the implicit sparse structure of the resulting Bayes net.
*
* This structure is examined even more closely in a JunctionTree, which
* additionally identifies cliques in the chordal Bayes net.
* \nosubgrouping
*/
struct ETNode
{
    int key; ///< key associated with root
    std::vector<int> factorsindex; ///< factors associated with root
    std::vector<int> childrenindex;
    int problemSize_;
    RealGaussianFactor* eliminate(GaussianBayesNetPointer* output,
                                  const int  Eliminatefunction,
                                  const std::vector<RealGaussianFactor*>& childrenFactors,const GaussianFactorGraph& gf);
};

class EliminationTree
{
public:
    std::vector<int> roots_;
    std::vector<ETNode*> nodes;
  //  std::vector<int> nodeskey;
    // std::vector<RealGaussianFactor> remainingFactors_;
    std::vector<int> remainingFactorsindex_;

    /// @name Standard Constructors
    /// @{

    /**
    * Build the elimination tree of a factor graph using pre-computed column structure.
    * @param factorGraph The factor graph for which to build the elimination tree
    * @param structure The set of factors involving each variable.  If this is not
    * precomputed, you can call the Create(const FactorGraph<DERIVEDFACTOR>&)
    * named constructor instead.
    * @return The elimination tree
    */
public:
    EliminationTree(const GaussianFactorGraph& factorGraph,
                           const VariableIndex& structure, const Ordering& order);

    /** Build the elimination tree of a factor graph.  Note that this has to compute the column
    * structure as a VariableIndex, so if you already have this precomputed, use the other
    * constructor instead.
    * @param factorGraph The factor graph for which to build the elimination tree
    */
    EliminationTree(const GaussianFactorGraph& factorGraph, const Ordering& order);

    /** Copy constructor - makes a deep copy of the tree structure, but only pointers to factors are
     *  copied, factors are not cloned. */
    EliminationTree(const EliminationTree& other)
    {

        roots_=other.roots_;
        nodes=other.nodes;
        //nodeskey=other.nodeskey;
        // Assign the remaining factors - these are pointers to factors in the original factor graph and
        // we do not clone them.
        remainingFactorsindex_ = other.remainingFactorsindex_;
    }

    /** Assignment operator - makes a deep copy of the tree structure, but only pointers to factors
     *  are copied, factors are not cloned. */
    EliminationTree& operator=(const EliminationTree& other);


    ~EliminationTree()
    {
       for(std::vector<ETNode*>::iterator bbi=nodes.begin();bbi!=nodes.end();bbi++)
       {
          delete *bbi;
          *bbi=NULL;
       }

    }

    /// @}

public:
    /// @name Standard Interface
    /// @{

    /** Eliminate the factors to a Bayes net and remaining factor graph
    * @param function The function to use to eliminate, see the namespace functions
    * in GaussianFactorGraph.h
    * @return The Bayes net and factor graph resulting from elimination
    */
    std::pair<GaussianBayesNetPointer*, GaussianFactorGraph*>
    eliminate(const int Eliminatefunction,const GaussianFactorGraph& gf) const;

    /// @}
    /// @name Testable
    /// @{

    /** Print the tree to cout */
    void print() const;

public:
    /// @name Advanced Interface
    /// @{
    /** Return the remaining factors that are not pulled into elimination */
    const std::vector<int> remainingFactorsindex() const
    {
        return remainingFactorsindex_;
    }

    /** Swap the data of this tree with another one, this operation is very fast. */
    void swap(EliminationTree& other);

protected:
    /// Protected default constructor
    EliminationTree() {}

};

struct ETEliminationDataChildrenFactorsPointer
{
    int parentindex_;
    int currentindex_;
    std::vector<RealGaussianFactor*> childFactors;
    ETEliminationDataChildrenFactorsPointer(
        int pIndex,int cindex,int nChildren):
        parentindex_(pIndex),currentindex_(cindex)
    {
        // childFactors=new std::vector<RealGaussianFactor>;
        childFactors.reserve(nChildren);
        // childFactors.resize(nChildren);
    }
};

struct ETTraversalNodePointer
{
    bool expanded;
    ETNode* treeNode;
    int parentindex;
// int currentpointerindex;

    //typename FastList<DATA>::iterator dataPointer;
    ETTraversalNodePointer(ETNode* _treeNode, int _parentindex) :
        expanded(false), treeNode(_treeNode), parentindex(_parentindex)
    {
    }
};

std::vector<RealGaussianFactor*>
inferenceEliminateTreePointer(GaussianBayesNetPointer* result, const EliminationTree& tree,
                              const int Eliminatefunction,const GaussianFactorGraph& gf);
//void ETDepthFirstForestClone(const EliminationTree& forest, ETNode* rootData);

void ETDepthFirstForestBNPointer(GaussianBayesNetPointer* BNresult,
                                 const int elimationtype,const EliminationTree& forest,
                                 ETEliminationDataChildrenFactorsPointer* rootData,const GaussianFactorGraph& gf);


#endif // ELIMINATIONTREEPOINTER_H_INCLUDED
