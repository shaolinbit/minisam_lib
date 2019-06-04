#ifndef LINEARCONTAINERFACTORPOSE_H_INCLUDED
#define LINEARCONTAINERFACTORPOSE_H_INCLUDED

#include "../nonlinear/NonlinearFactorGraph.h"
#include "../gmfconfig.h"

// Forward declarations
class JacobianFactor;
class HessianFactor;

/**
 * Dummy version of a generic linear factor to be injected into a nonlinear factor graph
 *
 * This factor does have the ability to perform relinearization under small-angle and
 * linearity assumptions if a linearization point is added.
 */
class  LinearContainerFactorPose : public NonlinearFactor
{
protected:

    RealGaussianFactor factor_;
    std::map<int,Eigen::VectorXd> linearizationPoint_;
#ifdef GMF_Using_Pose3
    std::map<int,Pose3> linearizationPose_;
#else
    std::map<int,Pose2> linearizationPose_;
#endif // GMF_Using_Pose3

public:
    /** Default constructor - necessary for serialization */
    LinearContainerFactorPose()
    {
        factortype=3;
    };
    ~LinearContainerFactorPose();

    /** direct copy constructor */
#ifdef GMF_Using_Pose3
    LinearContainerFactorPose(const RealGaussianFactor& factor, const std::map<int,Eigen::VectorXd>& linearizationPoint,
                              const  std::map<int,Pose3>& linearizationPose_):
        NonlinearFactor(factor.keys()), factor_(factor),
        linearizationPoint_(linearizationPoint),linearizationPose_(linearizationPose_)
    {
        factortype=3;
    }
    /** Primary constructor: store a linear factor with optional linearization point */
    LinearContainerFactorPose(const JacobianFactor& factor, const std::map<int,Eigen::VectorXd>& linearizationPoint,
                              const  std::map<int,Pose3>& linearizationPose): NonlinearFactor(factor.keys()),
        factor_(*(factor.clone()))
    {
        initializeLinearizationPoint(linearizationPoint);
        initializeLinearizationPose(linearizationPose);
        factortype=3;
    }

    /** Primary constructor: store a linear factor with optional linearization point */
    LinearContainerFactorPose(const HessianFactor& factor, const std::map<int,Eigen::VectorXd>& linearizationPoint,
                              const  std::map<int,Pose3>& linearizationPose):NonlinearFactor(factor.keys()),
        factor_(*(factor.clone()))
    {
        initializeLinearizationPoint(linearizationPoint);
        initializeLinearizationPose(linearizationPose);
        factortype=3;
    }
    /**
      * Calculate the nonlinear error for the factor, where the error is computed
      * by passing the delta between linearization point and c, where
      * delta = linearizationPoint_.localCoordinates(c), into the error function
      * of the stored linear factor.
      *
      * @return nonlinear error if linearizationPoint present, zero otherwise
      */
    double LCFPerrorP3v(const  std::map<int,Pose3>& linearizationPose_,const std::map<int,Eigen::VectorXd>& c) const;

    /** Extract the linearization point used in recalculating error */
    const std::map<int,Pose3>& linearizationPose() const
    {
        return linearizationPose_;
    }
#else
    LinearContainerFactorPose(const RealGaussianFactor& factor, const std::map<int,Eigen::VectorXd>& linearizationPoint,
                              const  std::map<int,Pose2>& linearizationPose_):
        NonlinearFactor(factor.keys()), factor_(factor),
        linearizationPoint_(linearizationPoint),
        linearizationPose_(linearizationPose_)
    {
        factortype=3;
    }
    /** Primary constructor: store a linear factor with optional linearization point */
    LinearContainerFactorPose(const JacobianFactor& factor, const std::map<int,Eigen::VectorXd>& linearizationPoint,
                              const  std::map<int,Pose2>& linearizationPose): NonlinearFactor(factor.keys()),
        factor_(factor)
    {
        initializeLinearizationPoint(linearizationPoint);
        initializeLinearizationPose(linearizationPose);
        factortype=3;
    }

    /** Primary constructor: store a linear factor with optional linearization point */
    LinearContainerFactorPose(const HessianFactor& factor, const std::map<int,Eigen::VectorXd>& linearizationPoint,
                              const  std::map<int,Pose2>& linearizationPose):NonlinearFactor(factor.keys()),
        factor_(factor)
    {
        initializeLinearizationPoint(linearizationPoint);
        initializeLinearizationPose(linearizationPose);
        factortype=3;
    }
    // Access
    /**
      * Calculate the nonlinear error for the factor, where the error is computed
      * by passing the delta between linearization point and c, where
      * delta = linearizationPoint_.localCoordinates(c), into the error function
      * of the stored linear factor.
      *
      * @return nonlinear error if linearizationPoint present, zero otherwise
      */
    double LCFPerrorP2v(const  std::map<int,Pose2>& Pose_,const std::map<int,Eigen::VectorXd>& c) const;
//virtual double LCFPerrorP2v(const std::map<int,Pose2> p3,const std::map<int,Eigen::VectorXd> v) const
    /** Extract the linearization point used in recalculating error */
    const std::map<int,Pose2>& linearizationPose() const
    {
        return linearizationPose_;
    }

#endif // GMF_Using_Pose3


    const RealGaussianFactor factor() const
    {
        return factor_;
    }



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
#ifdef GMF_Using_Pose3
    RealGaussianFactor linearizePV(const std::map<int,Pose3>& Pose_,
                                   const std::map<int,Eigen::VectorXd>& c) const;
    RealGaussianFactor* linearizePVPointer(const std::map<int,Pose3>& Pose_,
                                           const std::map<int,Eigen::VectorXd>& c,int factorization) const;
#else
    //RealGaussianFactor linearize(const std::map<int,Eigen::VectorXd>& c,const  std::map<int,Pose2>& Pose_) const;

    RealGaussianFactor linearizePV(const std::map<int,Pose2>& Pose_,
                                   const std::map<int,Eigen::VectorXd>& c) const;
    RealGaussianFactor* linearizePVPointer(const std::map<int,Pose2>& Pose_,
                                           const std::map<int,Eigen::VectorXd>& c,int factorization) const;
#endif // GMF_Using_Pose3

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
        LinearContainerFactorPose* nl=new LinearContainerFactorPose(factor_,linearizationPoint_,linearizationPose_);
        //  double bb=nl->LCFPerrorP2v(linearizationPose_,linearizationPoint_);
        return nl;
    }

    // casting syntactic sugar


    inline bool hasLinearizationPoint() const
    {
        return linearizationPoint_.size()>0;
    }

    inline bool hasLinearizationPose() const
    {
        return linearizationPose_.size()>0;
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
#ifdef GMF_Using_Pose3
    static NonlinearFactorGraph ConvertLinearGraph(const GaussianFactorGraph& linear_graph,

            const std::map<int,Eigen::VectorXd>& linearizationPoint,const  std::map<int,Pose3>& Pose_);

    void initializeLinearizationPose(const std::map<int,Pose3>& linearizationPose);
#else
    static NonlinearFactorGraph ConvertLinearGraph(const GaussianFactorGraph& linear_graph,

            const std::map<int,Eigen::VectorXd>& linearizationPoint,const  std::map<int,Pose2>& Pose_);
    void initializeLinearizationPose(const std::map<int,Pose2>& linearizationPose);
#endif // GMF_Using_Pose3
protected:
    void initializeLinearizationPoint(const std::map<int,Eigen::VectorXd>& linearizationPoint);


};


#endif // LINEARCONTAINERFACTORPOSE_H_INCLUDED
