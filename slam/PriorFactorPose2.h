#ifndef PRIORFACTORPOSE2_H_INCLUDED
#define PRIORFACTORPOSE2_H_INCLUDED



/**
 * A class for a soft prior on Pose2 Value
 * @addtogroup SLAM
 */

namespace minisam
{
class PriorFactorPose2: public NoiseModelFactor1
{
private:

    Pose2 prior_; /** The measurement */

public:

    /** default constructor - only use for serialization */
    PriorFactorPose2():NoiseModelFactor1(1) {}

    virtual ~PriorFactorPose2() {}

    /** Constructor */
    PriorFactorPose2(int key, const Pose2& prior,
                     GaussianNoiseModel* model) :
        NoiseModelFactor1(model, key,1), prior_(prior)
    {
    }

    /** Convenience constructor Pose3 takes a full covariance argument */
    PriorFactorPose2(int key, const Pose2& prior, const Eigen::MatrixXd& covariance) :
        NoiseModelFactor1(GaussianNoiseModel_Covariance(covariance), key,1), prior_(prior)
    {
    }

    virtual Eigen::VectorXd evaluateError(const Pose2& x) const
    {
        // Pose2 h=x.inverse()*prior_;
        Pose2 h=x.inverse()*prior_;
        Eigen::VectorXd v(3);
        v=-Pose2::ChartAtOrigin::Local(h);
        return v;
    }

    virtual Eigen::VectorXd evaluateError(const Pose2& x, Eigen::MatrixXd& H)const
    {
        H=Eigen::MatrixXd::Identity(3,3);
        //Pose2 h=x.inverse()*prior_;
        Pose2 h=x.inverse()*prior_;
        Eigen::VectorXd v(3);
        v=-Pose2::ChartAtOrigin::Local(h);
        return v;
    }


    const Pose2 & prior() const
    {
        return prior_;
    }

    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Pose2>& x,std::vector<Eigen::MatrixXd>& H) const
    {
        std::map<int,Pose2>::const_iterator itb=x.find(key());
        return evaluateError(itb->second,*(H.begin()));
    }

    //nonsense for virtual;
    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Eigen::VectorXd>& x)const
    {
        Eigen::VectorXd uw(3);
        uw.setZero();
        return uw;
    }
//nonsense for virtual;
    virtual Eigen::VectorXd evaluateError(const Eigen::VectorXd& x) const
    {

        Eigen::VectorXd uw(3);
        uw.setZero();
        return uw;

    }
    //nonsense for virtual;
    virtual Eigen::VectorXd evaluateError(const Eigen::VectorXd& x, Eigen::MatrixXd& H)const
    {
        Eigen::VectorXd uw(3);
        uw.setZero();
        return uw;
    }
    //nonsense for virtual;
    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Pose2>& x1,
                                            const std::map<int,Eigen::VectorXd>& x2,
                                            std::vector<Eigen::MatrixXd>& H) const
    {
        Eigen::VectorXd uw(3);
        uw.setZero();
        return uw;
    }

    virtual NonlinearFactor* clone()const
    {
        PriorFactorPose2* newfactor=new PriorFactorPose2(key(),prior(),noiseModel_);
        return newfactor;
    }

};
};
#endif // PRIORFACTORPOSE2_H_INCLUDED
