#ifndef NONLINEARFACTORGRAPH_H
#define NONLINEARFACTORGRAPH_H

/**
 * @file    NonlinearFactorGraph.h
 * @brief   Factor Graph Constsiting of non-linear factors
 */


#include "../nonlinear/NonlinearFactor.h"
#include "../inference/FactorGraph.h"
#include "../linear/HessianFactor.h"
#include "../linear/Scatter.h"
#include "../geometry/Pose3.h"
#include "../geometry/Pose2.h"
#include "../gmfconfig.h"
#include <functional>
namespace minisam
{

// Forward declarations
class GaussianFactorGraph;
class HessianFactor;
class Scatter;



/**
 * A non-linear factor graph is a graph of non-Gaussian, i.e. non-linear factors,
 * which derive from NonlinearFactor. The values structures are typically (in SAM) more general
 * than just vectors, e.g., Rot3 or Pose3, which are objects in non-linear manifolds.
 * Linearizing the non-linear factor graph creates a linear factor graph on the
 * tangent vector space at the linearization point. Because the tangent space is a true
 * vector space, the config type will be an VectorValues in that linearized factor graph.
 */
class  NonlinearFactorGraph: public FactorGraph<NoiseModelFactor>
{

public:
    /** Default constructor */
    NonlinearFactorGraph() {}

    /** Construct from iterator over factors */
    template<typename ITERATOR>
    NonlinearFactorGraph(ITERATOR firstFactor, ITERATOR lastFactor) : FactorGraph<NoiseModelFactor>(firstFactor, lastFactor) {}

    /** Construct from container of factors (shared_ptr or plain objects) */
    template<class CONTAINER>
    explicit NonlinearFactorGraph(const CONTAINER& factors) : FactorGraph<NoiseModelFactor>(factors) {}

    /** Implicit copy/downcast constructor to override explicit template container constructor */
    template<class DERIVEDFACTOR>
    NonlinearFactorGraph(const FactorGraph<DERIVEDFACTOR>& graph) : FactorGraph<NoiseModelFactor>(graph) {}

    NonlinearFactorGraph(const NonlinearFactorGraph& graph);

    void clearmemory();

    void clearall();


    /** unnormalized error, \f$ 0.5 \sum_i (h_i(X_i)-z)^2/\sigma^2 \f$ in the most common case */
    double error(const std::map<int,minimatrix*>& values) const;
    /// Linearize a nonlinear factor graph
    int linearize(const std::map<int,minimatrix*>& linearizationPoint,GaussianFactorGraph& lng,int factorization=0) const;

};

Scatter& scatterFromValues(const std::map<int,minivector>& values, std::vector<int>& ordering);
Scatter& scatterFromValues(const std::map<int,minivector>& values);
std::map<int,minivector> values_vectorXd_vector2d(const std::map<int,minivector>& valuesxd);

};
#endif // NonlinearFactorGraph_H
