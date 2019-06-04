#ifndef ELIMINATEABLEFACTORPOINTERGRAPH_H_INCLUDED
#define ELIMINATEABLEFACTORPOINTERGRAPH_H_INCLUDED

/**
 * @file    EliminateableFactorPointerGraph.h
 * @brief   Variable elimination algorithms for factor graphs
 * @author
 * @date
 */

#pragma once



#include "../inference/Ordering.h"
#include <utility>
#include "../linear/GaussianConditional.h"
#include "../linear/JacobianFactor.h"
#include "../linear/GaussianFactorGraph.h"
#include "../linear/GaussianBayesNetPointer.h"
#include "../inference/EliminationTree.h"
#include "../inference/VariableIndex.h"
#include "../inference/BayesTreePointer.h"
#include "../inference/JunctionTreePointer.h"

GaussianBayesNetPointer* eliminateSequential(const GaussianFactorGraph& gf,const Ordering& ordering,const int ElimateType);

GaussianBayesNetPointer* eliminateSequential(const GaussianFactorGraph& gf,const int ElimateType);

/** Do sequential elimination of the given \c variables in an ordering computed by COLAMD to
    *  produce a Bayes net and a remaining factor graph.  This computes the factorization \f$ p(X)
    *  = p(A|B) p(B) \f$, where \f$ A = \f$ \c variables, \f$ X \f$ is all the variables in the
    *  factor graph, and \f$ B = X\backslash A \f;*/
std::pair<GaussianBayesNetPointer*, GaussianFactorGraph*> eliminatePartialSequential(const std::vector<int>& variables,
        const GaussianFactorGraph& gf, int EliminationFunctionType) ;

std::pair<GaussianBayesNetPointer*, GaussianFactorGraph*> eliminatePartialSequential(const std::vector<int>& variables,
        const GaussianFactorGraph& gf, int EliminationFunctionType,const VariableIndex& variableIndex) ;


std::pair<GaussianBayesNetPointer*, GaussianFactorGraph*>  eliminatePartialSequential(
    const Ordering& ordering, const GaussianFactorGraph& gf,int EliminationFunctionType, const VariableIndex& variableIndex);

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
BayesTreePointer* eliminateMultifrontal(
    Ordering& ordering,
    const int Eliminatefunction,
    VariableIndex& variableIndex,const GaussianFactorGraph& gf);//0:QR,1:precholesky

BayesTreePointer* eliminateMultifrontal(
    Ordering& ordering,
    const int Eliminatefunction,
    const GaussianFactorGraph& gf);

BayesTreePointer* eliminateMultifrontal(
    const int Eliminatefunction,
    const GaussianFactorGraph& gf);

/** Compute the marginal of the requested variables and return the result as a Bayes net.
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
GaussianBayesNetPointer* marginalMultifrontalBayesNet(const GaussianFactorGraph& gf,
        Ordering& variables,
        const int Eliminatefunction);

/** Compute the marginal factor graph of the requested variables. */
GaussianFactorGraph* marginal(
    const std::vector<int>& variables,
    const  int EliminatefunctionType,
    const VariableIndex& variableIndex,const GaussianFactorGraph& gf);

/** Compute the marginal factor graph of the requested variables. */
GaussianFactorGraph* marginalgf(const GaussianFactorGraph& gf,
                                       const std::vector<int>& variables,
                                       const  int EliminatefunctionType);

std::pair<BayesTreePointer*, GaussianFactorGraph*>
eliminatePartialMultifrontalKey(const std::vector<int>& variables,
                                const int Eliminatefunction,const VariableIndex& variableIndex,const GaussianFactorGraph& gf) ;

std::pair<BayesTreePointer*, GaussianFactorGraph*>
eliminatePartialMultifrontalKey(
    const std::vector<int>& variables, const int Eliminatefunction,const GaussianFactorGraph& gf) ;
Ordering OrderingColamdConstrainedFirst(const GaussianFactorGraph& graph,
        const std::vector<int>& constrainFirst, bool forceOrder=false);

std::pair<BayesTreePointer*, GaussianFactorGraph*>
eliminatePartialMultifrontal(const Ordering& ordering,
                             const int Eliminatefunction,const VariableIndex& variableIndex,const GaussianFactorGraph& gf);

std::pair<BayesTreePointer*, GaussianFactorGraph*>
eliminatePartialMultifrontal(const Ordering& ordering,
                             const int Eliminatefunction,const GaussianFactorGraph& gf);



#endif // ELIMINATEABLEFACTORPOINTERGRAPH_H_INCLUDED
