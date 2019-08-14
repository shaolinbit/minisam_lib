#ifndef ISAM2DATA_H_INCLUDED
#define ISAM2DATA_H_INCLUDED
#include "../nonlinear/NonlinearFactorGraph.h"
namespace minisam
{
class ISAM2Data
{
public:
    /** The current linearization point */
    std::map<int,Eigen::VectorXd> theta_;
    std::map<int,Eigen::VectorXd> resulttheta_;

#ifdef GMF_Using_Pose3
    std::map<int,Pose3> thetaPose_;
    std::map<int,Pose3> resultPose_;
#else
    std::map<int,Pose2> thetaPose_;
    std::map<int,Pose2> resultPose_;
#endif // GMF_Using_Pose3

    std::map<int,Eigen::VectorXd> delta_;

    mutable std::map<int,Eigen::VectorXd> deltaNewton_; // Only used when using Dogleg - stores the Gauss-Newton update
    mutable std::map<int,Eigen::VectorXd> RgProd_; // Only used when using Dogleg - stores R*g and is updated incrementally
    NonlinearFactorGraph nonlinearFactors_;
    /** The current linear factors, which are only updated as needed */
    mutable GaussianFactorGraph linearFactors_;
    VariableIndex variableIndex_;
    mutable std::set<int> deltaReplacedMask_; // TODO: Make sure accessed in the right way
    /** Set of variables that are involved with linear factors from marginalized
     * variables and thus cannot have their linearization points changed. */
    std::set<int> fixedVariables_;
    ISAM2Data();
    ~ISAM2Data();

    void clearfactors();
    void clearpose();
};

};
#endif // ISAM2DATA_H_INCLUDED
