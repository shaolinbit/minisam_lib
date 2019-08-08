#ifndef UNARYFACTOR_H_INCLUDED
#define UNARYFACTOR_H_INCLUDED
#include "../nonlinear/NonlinearFactor.h"


namespace minisam
{


class UnaryFactor: public NoiseModelFactor1
{

    // The factor will hold a measurement consisting of an (X,Y) location
    // We could this with a Point2 but here we just use two doubles
public:
    double mx_, my_;

public:
    /// shorthand for a smart pointer to a factor
    //typedef boost::shared_ptr<UnaryFactor> shared_ptr;

    // The constructor requires the variable key, the (X, Y) measurement value, and the noise model
    UnaryFactor(int j, double x, double y, GaussianNoiseModel* model):
        NoiseModelFactor1(model, j,1), mx_(x), my_(y) {}

    virtual ~UnaryFactor() {}

    // Using the NoiseModelFactor1 base class there are two functions that must be overridden.
    // The first is the 'evaluateError' function. This function implements the desired measurement
    // function, returning a vector of errors when evaluated at the provided variable value. It
    // must also calculate the Jacobians for this measurement function, if requested.
    virtual Eigen::VectorXd evaluateError(const Pose2& x) const
    {
        Eigen::VectorXd br(2);
        br<< x.x() - mx_, x.y() - my_;
        return br;
        //return (Vector(2) << q.x() - mx_, q.y() - my_).finished();
    }
    virtual Eigen::VectorXd evaluateError(const Pose2& x, Eigen::MatrixXd& H)const
    {
        // The measurement function for a GPS-like measurement is simple:
        // error_x = pose.x - measurement.x
        // error_y = pose.y - measurement.y
        // Consequently, the Jacobians are:
        // [ derror_x/dx  derror_x/dy  derror_x/dtheta ] = [1 0 0]
        // [ derror_y/dx  derror_y/dy  derror_y/dtheta ] = [0 1 0]
        //Eigen::MatrixXd H(2,3);
        // H<< 1.0,0.0,0.0, 0.0,1.0,0.0;
        H= (Eigen::MatrixXd(2,3) << 1.0,0.0,0.0, 0.0,1.0,0.0).finished();
        Eigen::VectorXd br(2);
        br<< x.x() - mx_, x.y() - my_;
        return br;

    }

    /*

    virtual Eigen::VectorXd evaluateError(const Eigen::VectorXd x) const
    {
        Eigen::VectorXd x0;
        return x0;
    }
    virtual Eigen::VectorXd evaluateError(const Eigen::VectorXd x, Eigen::MatrixXd& H)const
    {
        Eigen::VectorXd x0;
        return x0;
    }
    */


    /*Eigen::VectorXd evaluateError(const Pose2& q, boost::optional<Matrix&> H = boost::none) const
    {
      // The measurement function for a GPS-like measurement is simple:
      // error_x = pose.x - measurement.x
      // error_y = pose.y - measurement.y
      // Consequently, the Jacobians are:
      // [ derror_x/dx  derror_x/dy  derror_x/dtheta ] = [1 0 0]
      // [ derror_y/dx  derror_y/dy  derror_y/dtheta ] = [0 1 0]
      if (H) (*H) = (Matrix(2,3) << 1.0,0.0,0.0, 0.0,1.0,0.0).finished();
      return (Vector(2) << q.x() - mx_, q.y() - my_).finished();
    }*/

    // The second is a 'clone' function that allows the factor to be copied. Under most
    // circumstances, the following code that employs the default copy constructor should
    // work fine.

    virtual NonlinearFactor* clone()const
    {
        UnaryFactor* newfactor=new UnaryFactor(key(),mx_,my_,noiseModel_);
        return newfactor;
    }
}; // UnaryFactor
};

#endif // UNARYFACTOR_H_INCLUDED
