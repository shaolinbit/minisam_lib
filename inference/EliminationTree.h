#ifndef ELIMINATIONTREE_H_INCLUDED
#define ELIMINATIONTREE_H_INCLUDED

/* ----------------------------------------------------------------------------

* GTSAM Copyright 2010, Georgia Tech Research Corporation,
* Atlanta, Georgia 30332-0415
* All Rights Reserved
* Authors: Frank Dellaert, et al. (see THANKS for the full author list)

* See LICENSE for the license information

* -------------------------------------------------------------------------- */

/**
* @file    EliminationTree.h
* @author  Frank Dellaert
* @author  Richard Roberts
* @date    Oct 13, 2010
*/
#include <utility>
#include "../linear/GaussianFactorGraph.h"
#include "../linear/RealGaussianFactor.h"


namespace minisam
{

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
};

class EliminationTree
{
public:
    std::vector<int> roots_;
    std::vector<ETNode*> nodes;
    std::vector<int> remainingFactorsindex_;

    /// @name Standard Constructors
    /// @{

    /**
    * Build the elimination tree of a factor graph using pre-computed column structure.
    * @param factorGraph The factor graph for which to build the elimination tree
    * @param structure The set of factors involving each variable.  If this is not
    * precomputed, you can call the Create(const GaussianFactorGraph&)
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
    EliminationTree(const EliminationTree& other);

    /** Assignment operator - makes a deep copy of the tree structure, but only pointers to factors
     *  are copied, factors are not cloned. */
    EliminationTree& operator=(const EliminationTree& other);


    ~EliminationTree();

    /// @}

public:
    /// @name Testable
    /// @{

    /** Print the tree to cout */
    void print() const;

    ///@}

public:
    /// @name Advanced Interface
    /// @{
    /** Return the remaining factors that are not pulled into elimination */
    const std::vector<int> remainingFactorsindex() const;

    /** Swap the data of this tree with another one, this operation is very fast. */
    void swap(EliminationTree& other);

    ///@}

protected:
    /// Protected default constructor
    EliminationTree() {}

};

};
#endif // ELIMINATIONTREE_H_INCLUDED
