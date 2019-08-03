#ifndef JUNCTIONTREE_H_INCLUDED
#define JUNCTIONTREE_H_INCLUDED

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file JunctionTree.h
 * @date Feb 4, 2010
 * @author Kai Ni
 * @author Frank Dellaert
 * @author Richard Roberts
 * @brief The junction tree
 */

#pragma once

#include "../inference/ClusterTree.h"
#include "../inference/EliminationTree.h"
#include "../symbolic/SymbolicConditional.h"
namespace minisam
{

// Forward declarations
class EliminationTree;
/**
 * A JunctionTree is a cluster tree, a set of variable clusters with factors, arranged in a tree,
 * with the additional property that it represents the clique tree associated with a Bayes Net.
 *
 * A junction tree is an intermediate data structure in multifrontal variable
 * elimination.  Each node is a cluster of factors, along with a clique of variables that are
 * eliminated all at once. In detail, every node k represents a clique (maximal fully connected
 * subset) of an associated chordal graph, such as a chordal Bayes net resulting from elimination.
 *
 * The difference with the BayesTree is that a JunctionTree stores factors, whereas a
 * BayesTree stores conditionals, that are the product of eliminating the factors in the
 * corresponding JunctionTree cliques.
 *
 * The tree structure and elimination method are exactly analagous to the EliminationTree,
 * except that in the JunctionTree, at each node multiple variables are eliminated at a time.
 *
 * \addtogroup Multifrontal
 * \nosubgrouping
 */
class JunctionTree : public EliminatableClusterTree
{
public:

    /// @name Standard Constructors
    /// @{

    /** Build the junction tree from an elimination tree. */
    JunctionTree(const EliminationTree& eliminationTree,const GaussianFactorGraph& gf);
    /// @}

private:

    // Private default constructor (used in static construction methods)
    JunctionTree() {}

};

// Pre-order visitor function
int ConstructorTraversalVisitorPre(
    const ETNode& node,
    int parentDataindex,int increment,ConstructorTraversalDataChildFactors& myData);
// Post-order visitor function

void ConstructorTraversalVisitorPostAlg2(
    const ETNode& ETreeNode,ConstructorTraversalDataChildFactors* currentData,
    ConstructorTraversalDataChildFactors* parentData,std::vector<Cluster*>& ctlist,
    const GaussianFactorGraph& gf) ;

struct JTTraversalNode
{
    bool expanded;
    ETNode* treeNode;
    int parentindex;
    Cluster* parentcluster_;
    std::list<ConstructorTraversalDataChildFactors>::iterator dataPointer;
    JTTraversalNode(ETNode* _treeNode, int _parentindex,Cluster* parentcluster) :
        expanded(false), treeNode(_treeNode), parentindex(_parentindex),parentcluster_(parentcluster)
    {
    }
};

void JTDepthFirstForest(const EliminationTree& forest, ConstructorTraversalDataChildFactors* rootData,
                               std::vector<Cluster*>& ctlist,const GaussianFactorGraph& gf);

std::pair<SymbolicConditional*, Factor*>
EliminateSymbolic(const FactorGraph<Factor>& factors, const Ordering& keys);
};


#endif // JUNCTIONTREE_H_INCLUDED
