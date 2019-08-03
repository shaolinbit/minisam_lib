#ifndef GAUSSIANFACTORGRAPH_H_INCLUDED
#define GAUSSIANFACTORGRAPH_H_INCLUDED

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    GaussianFactorGraph.h
 * @brief   Linear Factor Graph where all factors are Gaussians
 * @author  Kai Ni
 * @author  Christian Potthast
 * @author  Alireza Fathi
 * @author  Richard Roberts
 * @author  Frank Dellaert
 */

#include "../inference/FactorGraph.h"
#include "../linear/NoiseModel.h"
#include "../linear/RealGaussianFactor.h"
#include "../linear/JacobianFactor.h"
#include "../linear/HessianFactor.h"
#include "../inference/Ordering.h"
#include <list>
#include <Eigen/Core>
#include <map>
namespace minisam
{

class JacobianFactor;
class Ordering;
/* ************************************************************************* */
/**
 * A Linear Factor Graph is a factor graph where all factors are Gaussian, i.e.
 *   Factor == GaussianFactor
 *   VectorValues = A values structure of vectors
 * Most of the time, linear factor graphs arise by linearizing a non-linear factor graph.
 */
class  GaussianFactorGraph :
    public FactorGraph<RealGaussianFactor>
{
public:
    /** Default constructor */
    GaussianFactorGraph();
    /** Construct from container of factors (shared_ptr or plain objects) */
    template<class CONTAINER>
    explicit GaussianFactorGraph(const CONTAINER& factors) : FactorGraph<RealGaussianFactor>(factors) {}


    GaussianFactorGraph(const std::vector<RealGaussianFactor*> factors) : FactorGraph<RealGaussianFactor>(factors) {}

    /** Virtual destructor */
    ~GaussianFactorGraph();// {}

    void clear();
    void clearmemory();

    void clearall();

    /**
     * Return the set of variables involved in the factors (computes a set
     * union).
     */
    std::set<int> keys() const;

    /** unnormalized error */
    double error(const std::map<int,Eigen::VectorXd>& x) const;
    std::map<int,Eigen::VectorXd> optimizeGradientSearch() const;
   // std::map<int,Eigen::VectorXd> gradient(const std::map<int,Eigen::VectorXd>& x0) const;
    std::map<int,Eigen::VectorXd> gradientAtZero() const;
    std::list<Eigen::VectorXd> operator*(const std::map<int,Eigen::VectorXd>& x) const;
    std::map<int,Eigen::VectorXd> optimize(Ordering& ordering, const int Eliminatekind) const;
    std::map<int,Eigen::VectorXd> hessianDiagonal() const;

};

/**
 * Evaluates whether linear factors have any constrained noise models
 * @return true if any factor is constrained.
 */
bool hasConstraints(const std::vector<const RealGaussianFactor*>& factors);

std::map<int, std::vector<int>> VariableSlots(const GaussianFactorGraph& factorGraph);
std::map<int, std::vector<int>> VariableSlots(const std::vector<const RealGaussianFactor*>& factorGraph);

std::map<int, std::vector<int>> VariableSlots(std::vector<RealGaussianFactor>& factors1,
                             std::vector<RealGaussianFactor>& factors2);

JacobianFactor* convertToJacobianFactorPtr(RealGaussianFactor* gf);
const JacobianFactor* convertToJacobianFactorPtr(const RealGaussianFactor* gf);

};

#endif // GAUSSIANFACTORGRAPH_H_INCLUDED
