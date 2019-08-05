#ifndef NONLINEARFACTORGRAPH_H
#define NONLINEARFACTORGRAPH_H

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    NonlinearFactorGraph.h
 * @brief   Factor Graph Constsiting of non-linear factors
 * @author  Frank Dellaert
 * @author  Carlos Nieto
 * @author  Christian Potthast
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
class  NonlinearFactorGraph: public FactorGraph<NonlinearFactor>
{

public:
    /** Default constructor */
     NonlinearFactorGraph() {}

    /** Construct from iterator over factors */
    template<typename ITERATOR>
    NonlinearFactorGraph(ITERATOR firstFactor, ITERATOR lastFactor) : FactorGraph<NonlinearFactor>(firstFactor, lastFactor) {}

    /** Construct from container of factors (shared_ptr or plain objects) */
    template<class CONTAINER>
    explicit NonlinearFactorGraph(const CONTAINER& factors) : FactorGraph<NonlinearFactor>(factors) {}

    /** Implicit copy/downcast constructor to override explicit template container constructor */
    template<class DERIVEDFACTOR>
    NonlinearFactorGraph(const FactorGraph<DERIVEDFACTOR>& graph) : FactorGraph<NonlinearFactor>(graph) {}

    NonlinearFactorGraph(const NonlinearFactorGraph& graph);

    void clearmemory();

    void clearall();

    /** Write the graph in GraphViz format for visualization*/
    template<typename T>
    void saveGraph(std::ostream& stm, const std::map<int,T>& values,
                   const GraphvizFormatting& graphvizFormatting = GraphvizFormatting()) const;

    void saveGraph(std::ostream& stm, const std::map<int,Eigen::Vector2d>& valuesv2,const std::map<int,Pose2>& valuesp2,
                   const std::map<int,Eigen::Vector3d>& valuesv3,const std::map<int,Pose3>& valuesp3,
                   const GraphvizFormatting& graphvizFormatting = GraphvizFormatting()) const;

    /** unnormalized error, \f$ 0.5 \sum_i (h_i(X_i)-z)^2/\sigma^2 \f$ in the most common case */
    double error(const std::map<int,Eigen::VectorXd>& values) const;


#ifdef GMF_Using_Pose3
    /** unnormalized error, \f$ 0.5 \sum_i (h_i(X_i)-z)^2/\sigma^2 \f$ in the most common case */
    double error(const std::map<int,Eigen::VectorXd>& valuevector,const std::map<int,Pose3>& valuepose) const;
#else
    double error(const std::map<int,Eigen::VectorXd>& valuevector,const std::map<int,Pose2>& valuepose) const;
#endif // GMF_Using_Pose3

    /// Linearize a nonlinear factor graph
    int linearize(const std::map<int,Eigen::VectorXd>& linearizationPoint,GaussianFactorGraph& lng,int factorization=0) const;

#ifdef GMF_Using_Pose3
    int linearize(const std::map<int,Eigen::VectorXd>& linearizationPoint,const std::map<int,Pose3>& linearizationpose,
                         GaussianFactorGraph& lng,int factorization=0) const;
#else
    int linearize(const std::map<int,Eigen::VectorXd>& linearizationPoint,const std::map<int,Pose2>& linearizationpose,
                         GaussianFactorGraph& lng,int factorization) const;
#endif // GMF_Using_Pose3


};

Scatter& scatterFromValues(const std::map<int,Eigen::VectorXd>& values, std::vector<int>& ordering);
Scatter& scatterFromValues(const std::map<int,Eigen::VectorXd>& values);
std::map<int,Eigen::Vector2d> values_vectorXd_vector2d(const std::map<int,Eigen::VectorXd>& valuesxd);
void graphvizgetXY(Eigen::Vector2d* xygot,Eigen::Vector2d* valuev2,Pose2* valuep2,
                   Eigen::Vector3d* valuev3,Pose3* valuep3,  const GraphvizFormatting& graphvizFormatting);

};
#endif // NonlinearFactorGraph_H
