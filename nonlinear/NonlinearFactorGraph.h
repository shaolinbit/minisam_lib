#ifndef NONLINEARFACTORPOINTERGRAPH_H
#define NONLINEARFACTORPOINTERGRAPH_H


/**
 * @file    NonlinearFactorGraph.h
 * @brief   Factor Graph Constsiting of non-linear factors
 * @author
 */


#include "../nonlinear/NonlinearFactor.h"
#include "../inference/FactorPointerGraph.h"
#include "../linear/HessianFactor.h"
#include "../linear/Scatter.h"
#include "../geometry/Pose3.h"
#include "../geometry/Pose2.h"
#include "../gmfconfig.h"
#include <functional>


// Forward declarations
class Ordering;
//class GaussianFactorGraph;
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
class  NonlinearFactorGraph: public FactorPointerGraph<NonlinearFactor>
{

public:
    /** Default constructor */
    NonlinearFactorGraph() {}

    /** Construct from iterator over factors */
    template<typename ITERATOR>
    NonlinearFactorGraph(ITERATOR firstFactor, ITERATOR lastFactor) : FactorPointerGraph<NonlinearFactor>(firstFactor, lastFactor) {}

    /** Construct from container of factors (shared_ptr or plain objects) */
    template<class CONTAINER>
    explicit NonlinearFactorGraph(const CONTAINER& factors) : FactorPointerGraph<NonlinearFactor>(factors) {}

    /** Implicit copy/downcast constructor to override explicit template container constructor */
    template<class DERIVEDFACTOR>
    NonlinearFactorGraph(const FactorPointerGraph<DERIVEDFACTOR>& graph) : FactorPointerGraph<NonlinearFactor>(graph) {}

    NonlinearFactorGraph(const NonlinearFactorGraph& graph)
    {
        factors_.clear();
        for(auto& fg: graph)
        {
            factors_.push_back(fg);
        }
    }

    void clearmemory()
    {
       for(std::vector<NonlinearFactor*>::iterator ffi=this->begin();ffi!=this->end();ffi++)
       {
         delete *ffi;
         *ffi=NULL;
       }
    }
    void clearall()
    {
       this->clearmemory();
       this->clear();
    }

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

    /** Unnormalized probability. O(n) */
    double probPrime(const std::map<int,Eigen::VectorXd>& values) const;

    /**
     * Create a symbolic factor graph
     */
    FactorPointerGraph<Factor> symbolic();

    /**
     * Compute a fill-reducing ordering using COLAMD.
     */
    Ordering orderingCOLAMD() const;

    /**
     * Compute a fill-reducing ordering with constraints using CCOLAMD
     *
     * @param constraints is a map of Key->group, where 0 is unconstrained, and higher
     * group numbers are further back in the ordering. Only keys with nonzero group
     * indices need to appear in the constraints, unconstrained is assumed for all
     * other variables
     */
    Ordering orderingCOLAMDConstrained(const std::map<int, int>& constraints) const;

    /// Linearize a nonlinear factor graph
    //GaussianFactorGraph linearize(const std::map<int,Eigen::VectorXd>& linearizationPoint) const;

     int linearizePointer(const std::map<int,Eigen::VectorXd>& linearizationPoint,GaussianFactorGraph& lng,int factorization) const;

#ifdef GMF_Using_Pose3
    int linearizePointer(const std::map<int,Eigen::VectorXd>& linearizationPoint,const std::map<int,Pose3>& linearizationpose,
                         GaussianFactorGraph& lng,int factorization) const;
#else
    int linearizePointer(const std::map<int,Eigen::VectorXd>& linearizationPoint,const std::map<int,Pose2>& linearizationpose,
                         GaussianFactorGraph& lng,int factorization) const;
#endif // GMF_Using_Pose3

    /// typdef for dampen functions used below
    typedef std::function<void(const HessianFactor& hessianFactor)> Dampen;

    /**
     * Instead of producing a GaussianFactorGraph, pre-allocate and linearize directly
     * into a HessianFactor. Avoids the many mallocs and pointer indirection in constructing
     * a new graph, and hence useful in case a dense solve is appropriate for your problem.
     * An optional ordering can be given that still decides how the Hessian is laid out.
     * An optional lambda function can be used to apply damping on the filled Hessian.
     * No parallelism is exploited, because all the factors write in the same memory.
     */
    HessianFactor linearizeToHessianFactor(
        const std::map<int,Eigen::VectorXd>& values,
        const Dampen& dampen = nullptr) const;


    HessianFactor linearizeToHessianFactor(
        const std::map<int,Eigen::VectorXd>& values, Ordering& ordering,
        const Dampen& dampen = nullptr) const;


    /// Linearize and solve in one pass.
    /// Calls linearizeToHessianFactor, densely solves the normal equations, and updates the values.
    std::map<int,Eigen::VectorXd> updateCholesky(const std::map<int,Eigen::VectorXd>& values, Ordering& ordering,
            const Dampen& dampen = nullptr) const;

    std::map<int,Eigen::VectorXd> updateCholesky(const std::map<int,Eigen::VectorXd>& values,
            const Dampen& dampen = nullptr) const;

};

Scatter scatterFromValues(const std::map<int,Eigen::VectorXd>& values, Ordering& ordering);
Scatter scatterFromValues(const std::map<int,Eigen::VectorXd>& values);
std::map<int,Eigen::Vector2d> values_vectorXd_vector2d(const std::map<int,Eigen::VectorXd>& valuesxd);
void graphvizgetXY(Eigen::Vector2d* xygot,Eigen::Vector2d* valuev2,Pose2* valuep2,
                   Eigen::Vector3d* valuev3,Pose3* valuep3,  const GraphvizFormatting& graphvizFormatting);
#endif // NonlinearFactorGraph_H
