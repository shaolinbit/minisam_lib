#ifndef MESTIMATOR_H_INCLUDED
#define MESTIMATOR_H_INCLUDED
#include "../linear/NoiseModel.h"
namespace minisam
{
      class  Base_mEstimator {
      public:
        enum ReweightScheme { Scalar, Block };

      protected:
        /** the rows can be weighted independently according to the error
        * or uniformly with the norm of the right hand side */
        ReweightScheme reweight_;

      public:
        Base_mEstimator(const ReweightScheme reweight = Block):reweight_(reweight) {}
        virtual ~Base_mEstimator() {}

        /*
         * This method is responsible for returning the total penalty for a given amount of error.
         * For example, this method is responsible for implementing the quadratic function for an
         * L2 penalty, the absolute value function for an L1 penalty, etc.
         *
         */
        virtual double residual(double error) const { return 0; }

        /*
         * This method is responsible for returning the weight function for a given amount of error.
         * The weight function is related to the analytic derivative of the residual function. See
         *  http://research.microsoft.com/en-us/um/people/zhang/INRIA/Publis/Tutorial-Estim/node24.html
         * for details. This method is required when optimizing cost functions with robust penalties
         * using iteratively re-weighted least squares.
         */
        virtual double weight(double error) const{ return 0;}


        double sqrtWeight(double error) const {
          return std::sqrt(weight(error));
        }

        /** produce a weight vector according to an error vector and the implemented
        * robust function */
        minivector weight(const minivector &error) const;

        /** square root version of the weight function */
        minivector sqrtWeight(const minivector &error) const {
          //return weight(error).cwiseSqrt();
          minivector werror=weight(error);
          minivector_cwisesqrt(&werror);
          return werror;
        }

        /** reweight block matrices and a vector according to their weight implementation */
        void reweight(minivector &error) const;
        void reweight(std::vector<minimatrix> &A, minivector &error) const;
        void reweight(minimatrix &A, minivector &error) const;
        void reweight(minimatrix &A1, minimatrix &A2, minivector &error) const;
        void reweight(minimatrix &A1, minimatrix &A2, minimatrix &A3, minivector &error) const;

      };

      /// Null class is not robust so is a Gaussian ?
      class  Null_mEstimator : public Base_mEstimator {
      public:
        //typedef boost::shared_ptr<Null> shared_ptr;

        Null_mEstimator(const ReweightScheme reweight = Block) : Base_mEstimator(reweight) {}
        virtual ~Null_mEstimator() {}
        virtual double weight(double /*error*/) const { return 1.0; }
        static Null_mEstimator* Create() ;

      };

      /// Fair implements the "Fair" robust error model (Zhang97ivc)
      class  Fair_mEstimator : public Base_mEstimator {
      protected:
        double c_;

      public:
       // typedef boost::shared_ptr<Fair> shared_ptr;

        Fair_mEstimator(double c = 1.3998, const ReweightScheme reweight = Block);
        double weight(double error) const {
          return 1.0 / (1.0 + fabs(error) / c_);
        }
        static Fair_mEstimator* Create(double c, const ReweightScheme reweight = Block) ;

      };

      /// Huber implements the "Huber" robust error model (Zhang97ivc)
      class  Huber_mEstimator : public Base_mEstimator {
      protected:
        double k_;

      public:
       // typedef boost::shared_ptr<Huber> shared_ptr;

        Huber_mEstimator(double k = 1.345, const ReweightScheme reweight = Block);
        double weight(double error) const {
          return (error < k_) ? (1.0) : (k_ / fabs(error));
        }
        static Huber_mEstimator* Create(double k, const ReweightScheme reweight = Block) ;

      };

      /// Cauchy implements the "Cauchy" robust error model (Lee2013IROS).  Contributed by:
      ///   Dipl.-Inform. Jan Oberlaender (M.Sc.), FZI Research Center for
      ///   Information Technology, Karlsruhe, Germany.
      ///   oberlaender@fzi.de
      /// Thanks Jan!
      class  Cauchy_mEstimator : public Base_mEstimator {
      protected:
        double k_, ksquared_;

      public:
       // typedef boost::shared_ptr<Cauchy> shared_ptr;

        Cauchy_mEstimator(double k = 0.1, const ReweightScheme reweight = Block);
        double weight(double error) const {
          return ksquared_ / (ksquared_ + error*error);
        }
        static Cauchy_mEstimator* Create(double k, const ReweightScheme reweight = Block) ;

      };

      /// Tukey implements the "Tukey" robust error model (Zhang97ivc)
      class  Tukey_mEstimator : public Base_mEstimator {
      protected:
        double c_, csquared_;

      public:
        //typedef boost::shared_ptr<Tukey> shared_ptr;

        Tukey_mEstimator(double c = 4.6851, const ReweightScheme reweight = Block);
        double weight(double error) const {
          if (std::fabs(error) <= c_) {
            double xc2 = error*error/csquared_;
            return (1.0-xc2)*(1.0-xc2);
          }
          return 0.0;
        }
        static Tukey_mEstimator* Create(double k, const ReweightScheme reweight = Block) ;

      };

      /// Welsh implements the "Welsh" robust error model (Zhang97ivc)
      class  Welsh_mEstimator : public Base_mEstimator {
      protected:
        double c_, csquared_;

      public:
        //typedef boost::shared_ptr<Welsh> shared_ptr;

        Welsh_mEstimator(double c = 2.9846, const ReweightScheme reweight = Block);
        double weight(double error) const {
          double xc2 = (error*error)/csquared_;
          return std::exp(-xc2);
        }
        static Welsh_mEstimator* Create(double k, const ReweightScheme reweight = Block) ;

      };

      /// GemanMcClure implements the "Geman-McClure" robust error model
      /// (Zhang97ivc).
      ///
      /// Note that Geman-McClure weight function uses the parameter c == 1.0,
      /// but here it's allowed to use different values, so we actually have
      /// the generalized Geman-McClure from (Agarwal15phd).
      class  GemanMcClure_mEstimator : public Base_mEstimator {
      public:
        //typedef boost::shared_ptr<GemanMcClure> shared_ptr;

        GemanMcClure_mEstimator(double c = 1.0, const ReweightScheme reweight = Block);
        virtual ~GemanMcClure_mEstimator() {}
        virtual double weight(double error) const;
        static GemanMcClure_mEstimator* Create(double k, const ReweightScheme reweight = Block) ;

      protected:
        double c_;

      };

      /// DCS implements the Dynamic Covariance Scaling robust error model
      /// from the paper Robust Map Optimization (Agarwal13icra).
      ///
      /// Under the special condition of the parameter c == 1.0 and not
      /// forcing the output weight s <= 1.0, DCS is similar to Geman-McClure.
      class  DCS_mEstimator: public Base_mEstimator {
      public:
        //typedef boost::shared_ptr<DCS> shared_ptr;

        DCS_mEstimator(double c = 1.0, const ReweightScheme reweight = Block);
        virtual ~DCS_mEstimator() {}
        virtual double weight(double error) const;
        static DCS_mEstimator* Create(double k, const ReweightScheme reweight = Block) ;

      protected:
        double c_;

      };

      /// L2WithDeadZone implements a standard L2 penalty, but with a dead zone of width 2*k,
      /// centered at the origin. The resulting penalty within the dead zone is always zero, and
      /// grows quadratically outside the dead zone. In this sense, the L2WithDeadZone penalty is
      /// "robust to inliers", rather than being robust to outliers. This penalty can be used to
      /// create barrier functions in a general way.
      class  L2WithDeadZone_mEstimator : public Base_mEstimator {
      public:
          double k_;

      public:
          //typedef boost::shared_ptr<L2WithDeadZone> shared_ptr;

          L2WithDeadZone_mEstimator(double k, const ReweightScheme reweight = Block);
          double residual(double error) const {
            const double abs_error = fabs(error);
            return (abs_error < k_) ? 0.0 : 0.5*(k_-abs_error)*(k_-abs_error);
          }
          double weight(double error) const {
            // note that this code is slightly uglier than above, because there are three distinct
            // cases to handle (left of deadzone, deadzone, right of deadzone) instead of the two
            // cases (deadzone, non-deadzone) above.
            if (fabs(error) <= k_) return 0.0;
            else if (error > k_) return (-k_+error)/error;
            else return (k_+error)/error;
          }
          static L2WithDeadZone_mEstimator* Create(double k, const ReweightScheme reweight = Block);

      };
    /**
     *  Base class for robust error models
     *  The robust M-estimators above simply tell us how to re-weight the residual, and are
     *  isotropic kernels, in that they do not allow for correlated noise. They also have no way
     *  to scale the residual values, e.g., dividing by a single standard deviation.
     *  Hence, the actual robust noise model below does this scaling/whitening in sequence, by
     *  passing both a standard noise model and a robust estimator.
     *
     *  Taking as an example noise = Isotropic::Create(d, sigma),  we first divide the residuals
     *  uw = |Ax-b| by sigma by "whitening" the system (A,b), obtaining r = |Ax-b|/sigma, and
     *  then we pass the now whitened residual 'r' through the robust M-estimator.
     *  This is currently done by multiplying with sqrt(w), because the residuals will be squared
     *  again in error, yielding 0.5 \sum w(r)*r^2.
     *
     *  In other words, while sigma is expressed in the native residual units, a parameter like
     *  k in the Huber norm is expressed in whitened units, i.e., "nr of sigmas".
     */
    class  RobustNoiseModel : public GaussianNoiseModel {
   // public:
   //   typedef boost::shared_ptr<Robust> shared_ptr;

    public:
      //typedef mEstimator::Base RobustModel;
      //typedef noiseModel::Base NoiseModel;

      Base_mEstimator* robust_; ///< robust error function used
      GaussianNoiseModel* noise_;   ///< noise model used

    public:

      /// Default Constructor for serialization
      RobustNoiseModel() {};

      /// Constructor
      RobustNoiseModel( Base_mEstimator* robust, GaussianNoiseModel* noise)
      : GaussianNoiseModel(GaussianNoiseModel(noise->dim(),true)), robust_(robust), noise_(noise)
      {
         //invsigmas_=minivector(noise->dim(),1.0);
         minivector_resize(&invsigmas_,dim_);
         minivector_set_all(&invsigmas_,1.0);
         minivector_resize(&sigmas_,dim_);
         minivector_set_all(&sigmas_,1.0);
         //sigmas_=minivector(noise->dim(),1.0);
      }

      /// Destructor
      virtual ~RobustNoiseModel() {delete robust_;delete noise_;}

      /// Return the contained robust error function
      Base_mEstimator* robust() const { return robust_; }

      /// Return the contained noise model
      GaussianNoiseModel* noise()  { return noise_; }

      // TODO: functions below are dummy but necessary for the noiseModel::Base
      inline virtual minivector whiten(const minivector& v) const
      {
           minivector r(v);
           this->WhitenSystem(r);
            return r;
        }
      inline virtual minimatrix Whiten(const minimatrix& A) const
      {
          minivector b;
          minimatrix B(A);
           this->WhitenSystem(B,b);
           return B; }
      inline virtual minivector unwhiten(const minivector& /*v*/) const
      { throw std::invalid_argument("unwhiten is not currently supported for robust noise models."); }
      inline virtual double distance(const minivector& v) const
      {
         // return this->whiten(v).squaredNorm();
          double result;
          minivector b=this->whiten(v);
          miniblas_ddot(b,b,&result);
          return  result;
      }
      // TODO(mike): fold the use of the m-estimator residual(...) function into distance(...)
      inline virtual double distance_non_whitened(const minivector& v) const
      { return robust_->residual(norm2d(v)); }
      // TODO: these are really robust iterated re-weighting support functions
      virtual void WhitenSystem(minivector& b) const;
      virtual void WhitenSystem(std::vector<minimatrix>& A, minivector& b) const;
      virtual void WhitenSystem(minimatrix& A, minivector& b) const;
      virtual void WhitenSystem(minimatrix& A1, minimatrix& A2, minivector& b) const;
      virtual void WhitenSystem(minimatrix& A1, minimatrix& A2, minimatrix& A3, minivector& b) const;

      static RobustNoiseModel* Create(
        Base_mEstimator* robust, GaussianNoiseModel* noise);

    };

};


#endif // MESTIMATOR_H_INCLUDED
