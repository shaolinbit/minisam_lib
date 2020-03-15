#ifndef GAUSSIANFACTORGRAPH_H_INCLUDED
#define GAUSSIANFACTORGRAPH_H_INCLUDED

/**
 * @file    GaussianFactorGraph.h
 * @brief   Linear Factor Graph where all factors are Gaussians
 */

#include "../inference/FactorGraph.h"
#include "../linear/NoiseModel.h"
#include "../linear/RealGaussianFactor.h"
#include "../linear/JacobianFactor.h"
#include "../linear/HessianFactor.h"
#include <list>
#include <map>
namespace minisam
{

class JacobianFactor;
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
    double error(const std::map<int,minivector>& x) const;
    std::map<int,minivector> optimizeGradientSearch() const;
    std::map<int,minivector> gradientAtZero() const;
    std::list<minivector> operator*(const std::map<int,minivector>& x) const;
    std::map<int,minivector> optimize(std::vector<int>& ordering,
    const Factorization Eliminatekind=CHOLESKY) const;
    std::map<int,minivector> hessianDiagonal() const;

};

/**
 * Evaluates whether linear factors have any constrained noise models
 * @return true if any factor is constrained.
 */
bool hasConstraints(const std::vector<const RealGaussianFactor*>& factors);

std::map<int, std::vector<int>> VariableSlots(const GaussianFactorGraph& factorGraph);
std::map<int, std::vector<int>> VariableSlots(const std::vector<const RealGaussianFactor*>& factorGraph);

JacobianFactor* convertToJacobianFactorPtr(RealGaussianFactor* gf);
const JacobianFactor* convertToJacobianFactorPtr(const RealGaussianFactor* gf);

};

#endif // GAUSSIANFACTORGRAPH_H_INCLUDED
