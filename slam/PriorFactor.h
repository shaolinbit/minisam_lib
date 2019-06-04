
/**
 *  @file  PriorFactor.h
 *  @author
 **/
#pragma once

#include "../nonlinear/NonlinearFactor.h"
#include "../linear/NoiseModel.h"


/**
 * A class for a soft prior on any Value type
 * @addtogroup SLAM
 */
//template<class VALUE>
class PriorFactor: public NoiseModelFactor1
{

public:
    //typedef VALUE T;

private:
    Eigen::VectorXd prior_; /** The measurement */
public:
    /** default constructor - only use for serialization */
    PriorFactor() {}

    virtual ~PriorFactor() {}

    /** Constructor */
    PriorFactor(int key, const Eigen::VectorXd& prior, SharedNoiseModel* model) :
        NoiseModelFactor1(model, key), prior_(prior)
    {
    }

    /** Convenience constructor that takes a full covariance argument */
    PriorFactor(int key, const Eigen::VectorXd& prior, const Eigen::MatrixXd& covariance) :
        NoiseModelFactor1(newGaussianNoiseModelCovariance(covariance), key), prior_(prior)
    {}


    /** implement functions needed for Testable */

    /** print
    virtual void print(const std::string& s, const KeyFormatter& keyFormatter = DefaultKeyFormatter) const {
      std::cout << s << "PriorFactor on " << keyFormatter(this->key()) << "\n";
      traits<T>::Print(prior_, "  prior mean: ");
      this->noiseModel_->print("  noise model: ");
    }*/

    /** equals
    virtual bool equals(const NonlinearFactor& expected, double tol=1e-9) const {
      const This* e = dynamic_cast<const This*> (&expected);
      return e != NULL && Base::equals(*e, tol) && traits<T>::Equals(prior_, e->prior_, tol);
    } */

    /** implement functions needed to derive from Factor */

    /** vector of errors
    Vector evaluateError(const T& x, boost::optional<Matrix&> H = boost::none) const
     {
      if (H) (*H) = Matrix::Identity(traits<T>::GetDimension(x),traits<T>::GetDimension(x));
      // manifold equivalent of z-x -> Local(x,z)
      // TODO(ASL) Add Jacobians.
      return -traits<T>::Local(x, prior_);
    }*/

    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Eigen::VectorXd>& x)const
    {
        //std::map<int,Eigen::VectorXd>::const_iterator itb=x.begin();
        std::map<int,Eigen::VectorXd>::const_iterator itb=x.find(key());
        return -(prior_-itb->second);//evaluateError(x);(prior_-x);
    }
    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Eigen::VectorXd>& x,std::vector<Eigen::MatrixXd>& H) const
    {
        //std::map<int,Eigen::VectorXd>::const_iterator itb=x.begin();
        std::map<int,Eigen::VectorXd>::const_iterator itb=x.find(key());
        // return (prior_-itb->second);//evaluateError(x);(prior_-x);
        return evaluateError(itb->second,*(H.begin()));
    }

    Eigen::VectorXd evaluateError(const Eigen::VectorXd& x) const
    {
        // if (H) (*H) = Matrix::Identity(traits<T>::GetDimension(x),traits<T>::GetDimension(x));
        // manifold equivalent of z-x -> Local(x,z)
        // TODO(ASL) Add Jacobians.
        return -(prior_-x);
    }
    Eigen::VectorXd evaluateError(const Eigen::VectorXd& x, Eigen::MatrixXd& H) const
    {
        // H = Eigen::.rows()),x.cols());
        H=Eigen::MatrixXd::Identity(x.rows(),x.rows());

        //return -traits<T>::Local(x, prior_);
        //return -return x.localCoordinates(prior_);

        //Class h = x.between(prior_); // derivatives inlined below
        //TangentVector v = Class::ChartAtOrigin::Local(h);
        Eigen::VectorXd h(x.rows());
        h=prior_-x;
        // Eigen::VectorXd v(x.rows());

        return -h;
        // manifold equivalent of z-x -> Local(x,z)
        // TODO(ASL) Add Jacobians.
        //return (prior_-x);
    }




    const Eigen::VectorXd prior() const
    {
        return prior_;
    }


    virtual Eigen::VectorXd unwhitenedError(const std::map<int,Pose3>& x1,
                                            const std::map<int,Eigen::VectorXd>& x2,
                                            std::vector<Eigen::MatrixXd>& H) const
    {
        Eigen::VectorXd uw(6);
        uw.setZero();
        return uw;
    }

#ifdef GMF_Using_Pose3
    //nonsense for virtual;
    virtual Eigen::VectorXd evaluateError(const Pose3 x) const
    {
        Eigen::VectorXd pp(6);
        pp.setZero();
        return pp;
    }
    //nonsense for virtual;
    virtual Eigen::VectorXd evaluateError(const Pose3 x, Eigen::MatrixXd& H) const
    {
        Eigen::VectorXd pp(6);
        pp.setZero();
        return pp;

    }
#else
    virtual Eigen::VectorXd evaluateError(const Pose2& x) const
    {
        Eigen::VectorXd xb;
        return xb;
    }
    virtual Eigen::VectorXd evaluateError(const Pose2& x, Eigen::MatrixXd& H) const
    {

        Eigen::VectorXd xb;
        return xb;

    }
#endif // GMF_Using_Pose3
    virtual NonlinearFactor* clone()const
    {
        PriorFactor* newfactor=new PriorFactor(key(),prior(),noiseModel_);
        return newfactor;
    }
};


