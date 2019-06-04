#ifndef ISAM2JUNCTIONTREEPOINTER_H_INCLUDED
#define ISAM2JUNCTIONTREEPOINTER_H_INCLUDED

/**
 * @file JunctionTree.h
 * @date
 * @author
 * @brief The junction tree
 */

//#pragma once

#include "../inference/ClusterTree.h"
#include "../inference/EliminationTree.h"
#include "../inference/Conditional.h"
#include "../symbolic/SymbolicConditional.h"


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
// template<class BAYESTREE, class GRAPH>
class ISAM2JunctionTreePointer : public EliminatableClusterTree
{
public:

    /// @name Standard Constructors
    /// @{

    /** Build the junction tree from an elimination tree. */
    //template<class ETREE>
    /*static ISAM2JunctionTreePointer FromEliminationTree(const EliminationTree& eliminationTree,const GaussianFactorGraph& gf)
    {
        return ISAM2JunctionTreePointer(eliminationTree,gf);
    }*/

    /** Build the junction tree from an elimination tree. */
    // template<class ETREE_BAYESNET, class ETREE_GRAPH>
    ISAM2JunctionTreePointer(const EliminationTree& eliminationTree,const GaussianFactorGraph& gf);

    /// @}

private:

    // Private default constructor (used in static construction methods)
    ISAM2JunctionTreePointer() {}

};

std::pair<ISAM2*, GaussianFactorGraph*>
eliminatePartialMultifrontalISAM2Pointer(const Ordering& ordering,
        const int Eliminatefunction,const VariableIndex& variableIndex,const GaussianFactorGraph& gf);

std::pair<ISAM2*, GaussianFactorGraph*>
eliminatePartialMultifrontalISAM2Pointer(const Ordering& ordering,
        const int Eliminatefunction,const GaussianFactorGraph& gf);

/** Compute the marginal factor graph of the requested variables. */
GaussianFactorGraph* marginalISAM2Pointer(
    const std::vector<int>& variables,
    const  int EliminatefunctionType,
    const VariableIndex& variableIndex,const GaussianFactorGraph& gf);

/** Compute the marginal factor graph of the requested variables. */
GaussianFactorGraph* marginalgfISAM2Pointer(const GaussianFactorGraph& gf,
        const std::vector<int>& variables,
        const  int EliminatefunctionType);


/** Do multifrontal elimination of all variables to produce a Bayes tree.  If an ordering is not
   *  provided, the ordering will be computed using either COLAMD or METIS, dependeing on
   *  the parameter orderingType (Ordering::COLAMD or Ordering::METIS)
   *
   *  <b> Example - Full Cholesky elimination in COLAMD order: </b>
   *  \code
   *  boost::shared_ptr<GaussianBayesTree> result = graph.eliminateMultifrontal(EliminateCholesky);
   *  \endcode
   *
   *  <b> Example - Full QR elimination in specified order:
   *  \code
   *  boost::shared_ptr<GaussianBayesTree> result = graph.eliminateMultifrontal(EliminateQR, myOrdering);
   *  \endcode
   *
   *  <b> Example - Reusing an existing VariableIndex to improve performance, and using COLAMD ordering: </b>
   *  \code
   *  VariableIndex varIndex(graph); // Build variable index
   *  Data data = otherFunctionUsingVariableIndex(graph, varIndex); // Other code that uses variable index
   *  boost::shared_ptr<GaussianBayesTree> result = graph.eliminateMultifrontal(EliminateQR, boost::none, varIndex);
   *  \endcode
   *  */
ISAM2* eliminateMultifrontalISAM2Pointer(
    const Ordering& ordering,
    const int Eliminatefunction,
    const VariableIndex& variableIndex,const GaussianFactorGraph& gf);//0:QR,1:precholesky

ISAM2* eliminateMultifrontalISAM2Pointer(
    const Ordering& ordering,
    const int Eliminatefunction,
    const GaussianFactorGraph& gf);


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
ISAM2* marginalMultifrontalISAM2Pointer(
    const std::vector<int>& variables,
    const int Defaultkind,const GaussianFactorGraph& gf);

//std::pair<SymbolicConditional, Factor>
//    EliminateSymbolic(const FactorGraph<Factor>& factors, const Ordering& keys);


#endif // ISAM2JUNCTIONTREEPOINTER_H_INCLUDED
