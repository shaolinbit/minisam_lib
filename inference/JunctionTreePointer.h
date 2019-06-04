#ifndef JUNCTIONTREEPOINTER_H_INCLUDED
#define JUNCTIONTREEPOINTER_H_INCLUDED



/**
 * @file JunctionTree.h
 * @date
 * @author
 * @brief The junction tree

 */

#pragma once

#include "../inference/ClusterTree.h"
#include "../inference/EliminationTree.h"
#include "../inference/Conditional.h"
#include "../symbolic/SymbolicConditional.h"

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
class JunctionTreePointer : public EliminatableClusterTree
{
public:

    /// @name Standard Constructors
    /// @{

    /** Build the junction tree from an elimination tree. */
    /*static JunctionTree FromEliminationTree(const EliminationTree& eliminationTree,const GaussianFactorGraph& gf)
    {
        return JunctionTree(eliminationTree,gf);
    }*/

    /** Build the junction tree from an elimination tree. */
    JunctionTreePointer(const EliminationTree& eliminationTree,const GaussianFactorGraph& gf);
    /// @}

private:

    // Private default constructor (used in static construction methods)
    JunctionTreePointer() {}

};

// Pre-order visitor function
static int ConstructorTraversalVisitorPrePointer(
    const ETNode& node,
    int parentDataindex,int increment,ConstructorTraversalDataChildFactorsPointer& myData);
// Post-order visitor function

void ConstructorTraversalVisitorPostAlg2(
    const ETNode& ETreeNode,ConstructorTraversalDataChildFactorsPointer* currentData,
    ConstructorTraversalDataChildFactorsPointer* parentData,std::vector<Cluster*>& ctlist,
    const GaussianFactorGraph& gf) ;

struct JTTraversalNodePointer
{
    bool expanded;
    ETNode* treeNode;
    int parentindex;
    Cluster* parentcluster_;
    std::list<ConstructorTraversalDataChildFactorsPointer>::iterator dataPointer;
    JTTraversalNodePointer(ETNode* _treeNode, int _parentindex,Cluster* parentcluster) :
        expanded(false), treeNode(_treeNode), parentindex(_parentindex),parentcluster_(parentcluster)
    {
    }
};


/** Compute the marginal of the requested variables and return the result as a Bayes tree.
   *  @param variables Determines the variables whose marginal to compute, if provided as an
   *         Ordering they will be ordered in the returned BayesNet as specified, and if provided
   *         as a vector<Key> they will be ordered using constrained COLAMD.
   *  @param marginalizedVariableOrdering Optional ordering for the variables being marginalized
   *         out, i.e. all variables not in \c variables.  If this is boost::none, the ordering
   *         will be computed with COLAMD.
   *  @param function Optional dense elimination function, if not provided the default will be
   *         used.
   *  @param variableIndex Optional pre-computed VariableIndex for the factor graph, if not
   *         provided one will be computed. */
BayesTreePointer* marginalMultifrontalBayesTreePointer(
    const std::vector<int>& variables,
    const int Defaultkind,const GaussianFactorGraph& gf);

void JTDepthFirstForestPointer(const EliminationTree& forest, ConstructorTraversalDataChildFactorsPointer* rootData,
                               std::vector<Cluster*>& ctlist,const GaussianFactorGraph& gf);

std::pair<SymbolicConditional*, Factor*>
EliminateSymbolicPointer(const FactorPointerGraph<Factor>& factors, const Ordering& keys);



#endif // JUNCTIONTREEPOINTER_H_INCLUDED
