#ifndef LINEARCONTAINERFACTOR_H
#define LINEARCONTAINERFACTOR_H

/**
 * @file LinearContainerFactor.h
 *
 * @brief Wrap Jacobian and Hessian linear factors to allow simple injection into a nonlinear graph
 *
 * @date
 * @author
 */

#include "../nonlinear/NonlinearFactorGraph.h"


// Forward declarations
class JacobianFactor;
class HessianFactor;

/**
 * Dummy version of a generic linear factor to be injected into a nonlinear factor graph
 *
 * This factor does have the ability to perform relinearization under small-angle and
 * linearity assumptions if a linearization point is added.
 */
class  LinearContainerFactor : public NonlinearFactor
{
protected:

    RealGaussianFactor factor_;
    std::map<int,Eigen::VectorXd> linearizationPoint_;

public:
    /** Default constructor - necessary for serialization */
    LinearContainerFactor();
    ~LinearContainerFactor();

    /** direct copy constructor */
    LinearContainerFactor(const RealGaussianFactor& factor, const std::map<int,Eigen::VectorXd>& linearizationPoint);

    /** Primary constructor: store a linear factor with optional linearization point */
    LinearContainerFactor(const JacobianFactor& factor, const std::map<int,Eigen::VectorXd>& linearizationPoint);

    /** Primary constructor: store a linear factor with optional linearization point */
    LinearContainerFactor(const HessianFactor& factor, const std::map<int,Eigen::VectorXd>& linearizationPoint);
    // Access

    const RealGaussianFactor factor() const
    {
        return factor_;
    }


    /**
     * Calculate the nonlinear error for the factor, where the error is computed
     * by passing the delta between linearization point and c, where
     * delta = linearizationPoint_.localCoordinates(c), into the error function
     * of the stored linear factor.
     *
     * @return nonlinear error if linearizationPoint present, zero otherwise
     */
    double error(const std::map<int,Eigen::VectorXd>& c) const;

    /** get the dimension of the factor: rows of linear factor */
    int dim() const;

    /** Extract the linearization point used in recalculating error */
    const std::map<int,Eigen::VectorXd>& linearizationPoint() const
    {
        return linearizationPoint_;
    }

    /**
     * Linearize to a GaussianFactor, with method depending on the presence of a linearizationPoint
     *  - With no linearization point, returns a cloned version of the stored linear factor.
     *  - With a linearization point provided, returns a relinearized version of
     *  the linearized factor.
     *
     * The relinearization approach used computes a linear delta between the original linearization
     * point and the new values c, where delta = linearizationPoint_.localCoordinates(c), and
     * substitutes this change into the system.  This substitution is only really valid for
     * linear variable manifolds, and for any variables based on a non-commutative
     * manifold (such as Pose2, Pose3), the relinearized version will be effective
     * for only small angles.
     *
     * TODO: better approximation of relinearization
     * TODO: switchable modes for approximation technique
     */
    RealGaussianFactor linearize(const std::map<int,Eigen::VectorXd>& c) const;

    /**
     * Creates an anti-factor directly
     */
    RealGaussianFactor negateToGaussian() const;

    /**
     * Creates the equivalent anti-factor as another LinearContainerFactor.
     */
    NonlinearFactor negateToNonlinear() const;

    /**
     * Creates a shared_ptr clone of the factor - needs to be specialized to allow
     * for subclasses
     *
     * Clones the underlying linear factor
     */
    NonlinearFactor* clone() const
    {
        LinearContainerFactor* nl=new LinearContainerFactor(factor_,linearizationPoint_);
        return nl;
    }

    // casting syntactic sugar


    inline bool hasLinearizationPoint() const
    {
        return linearizationPoint_.size()>0;
    }

    /**
     * Simple checks whether this is a Jacobian or Hessian factor
     */
    bool isJacobian() const;
    bool isHessian() const;

    /** Casts to JacobianFactor */
    JacobianFactor toJacobian() const;

    /** Casts to HessianFactor */
    HessianFactor toHessian() const;

    /**
     * Utility function for converting linear graphs to nonlinear graphs
     * consisting of LinearContainerFactors.
     */
    static NonlinearFactorGraph ConvertLinearGraph(const GaussianFactorGraph& linear_graph,

            const std::map<int,Eigen::VectorXd>& linearizationPoint);


protected:
    void initializeLinearizationPoint(const std::map<int,Eigen::VectorXd>& linearizationPoint);

}; // \class LinearContainerFactor

#endif // LINEARCONTAINERFACTOR_H
