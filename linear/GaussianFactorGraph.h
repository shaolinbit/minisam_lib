#ifndef GAUSSIANFACTORPOINTERGRAPH_H_INCLUDED
#define GAUSSIANFACTORPOINTERGRAPH_H_INCLUDED

/**
 * @file    GaussianFactorGraph.h
 * @brief   Linear Factor Graph where all factors are Gaussians
 * @author
 */

#include "../inference/FactorPointerGraph.h"
#include "../linear/NoiseModel.h"
#include "../linear/RealGaussianFactor.h"
#include "../linear/JacobianFactor.h"
#include "../linear/HessianFactor.h"
#include "../inference/Ordering.h"
#include <list>
#include <Eigen/Core>
#include <map>

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
    public FactorPointerGraph<RealGaussianFactor>
{
public:
    /** Default constructor */
    GaussianFactorGraph();
    /** Construct from container of factors (shared_ptr or plain objects) */
    template<class CONTAINER>
    explicit GaussianFactorGraph(const CONTAINER& factors) : FactorPointerGraph<RealGaussianFactor>(factors) {}


    GaussianFactorGraph(const std::vector<RealGaussianFactor*> factors) : FactorPointerGraph<RealGaussianFactor>(factors) {}

    /** Virtual destructor */
    ~GaussianFactorGraph();// {}

    void clear()
    {
        this->factors_.clear();
    }
    void clearmemory()
    {
       for(std::vector<RealGaussianFactor*>::iterator ffi=this->begin();ffi!=this->end();ffi++)
       {
         delete *ffi;
         *ffi=NULL;
       }
    }

    void clearall()
    {
     clearmemory();
     clear();
    }


    /** Add a factor by value - makes a copy */
    void add(RealGaussianFactor* factor)
    {
        this->push_back((factor->clone()));
        //  push_back(factor.clone());
    }

    /** Add a factor by pointer - stores pointer without copying the factor*/
    //void add(JacobianFactor* factor);// { push_back(factor); }

    void add(const Eigen::VectorXd& b);// {
    //   add(JacobianFactor(b)); }

    void add(int key1, const Eigen::MatrixXd& A1,
             const Eigen::VectorXd& b,  DiagonalNoiseModel* model =new DiagonalNoiseModel()) ;//{
    //add(JacobianFactor(key1,A1,b,model)); }

    void add(int key1, const Eigen::MatrixXd& A1,
             int key2, const Eigen::MatrixXd& A2,
             const Eigen::VectorXd& b, DiagonalNoiseModel* model =new DiagonalNoiseModel());// {
    // add(JacobianFactor(key1,A1,key2,A2,b,model)); }

    void add(int key1, const Eigen::MatrixXd& A1,
             int key2, const Eigen::MatrixXd& A2,
             int key3, const Eigen::MatrixXd& A3,
             const Eigen::VectorXd& b, DiagonalNoiseModel* model =new  DiagonalNoiseModel());// {
    //add(JacobianFactor(key1,A1,key2,A2,key3,A3,b,model)); }
    template<class TERMS>
    void add(const TERMS& terms, const Eigen::VectorXd &b,  DiagonalNoiseModel* model =new DiagonalNoiseModel());

    // {
    //   add(JacobianFactor(terms,b,model)); }

    /**
     * Return the set of variables involved in the factors (computes a set
     * union).
     */
    // typedef KeySet Keys;
    std::set<int> keys() const;

    /* return a map of (Key, dimension) */
    std::map<int, int> getKeyDimMap() const;
    //void setRealGaussianFactor(int i,RealGaussianFactor* gf);

    /** unnormalized error */
    double error(const std::map<int,Eigen::VectorXd>& x) const
    {
        double total_error = 0.;
        for(const RealGaussianFactor* factor: *this)
        {
            if(factor->size()!=0)
                total_error += factor->error(x);
        }
        return total_error;
    }

    /** Unnormalized probability. O(n) */
    double probPrime(const std::map<int,Eigen::VectorXd>& c) const
    {
        return exp(-0.5 * error(c));
    }

    /**
     * Clone() performs a deep-copy of the graph, including all of the factors.
     * Cloning preserves null factors so indices for the original graph are still
     * valid for the cloned graph.
     */
    GaussianFactorGraph clone() const;

    /**
     * CloneToPtr() performs a simple assignment to a new graph and returns it.
     * There is no preservation of null factors!
     */
//    virtual GaussianFactorGraph& cloneToPtr() const;

    /**
     * Returns the negation of all factors in this graph - corresponds to antifactors.
     * Will convert all factors to HessianFactors due to negation of information.
     * Cloning preserves null factors so indices for the original graph are still
     * valid for the cloned graph.
     */
    GaussianFactorGraph negate() const;

    ///@name Linear Algebra
    ///@{

    /**
     * Return vector of i, j, and s to generate an m-by-n sparse Jacobian matrix,
     * where i(k) and j(k) are the base 0 row and column indices, s(k) a double.
     * The standard deviations are baked into A and b
     */
    //std::vector<std::tuple<int, int, double> > sparseJacobian() const;

    /**
     * Matrix version of sparseJacobian: generates a 3*m matrix with [i,j,s] entries
     * such that S(i(k),j(k)) = s(k), which can be given to MATLAB's sparse.
     * The standard deviations are baked into A and b
     */
    // Eigen::MatrixXd sparseJacobian_() const;

    /**
     * Return a dense \f$ [ \;A\;b\; ] \in \mathbb{R}^{m \times n+1} \f$
     * Jacobian Eigen::MatrixXd, augmented with b with the noise models baked
     * into A and b.  The negative log-likelihood is
     * \f$ \frac{1}{2} \Vert Ax-b \Vert^2 \f$.  See also
     * GaussianFactorGraph::jacobian and GaussianFactorGraph::sparseJacobian.
     */
    Eigen::MatrixXd augmentedJacobian(const Ordering& optionalOrdering) const;
    Eigen::MatrixXd augmentedJacobian() const;

    /**
     * Return the dense Jacobian \f$ A \f$ and right-hand-side \f$ b \f$,
     * with the noise models baked into A and b. The negative log-likelihood
     * is \f$ \frac{1}{2} \Vert Ax-b \Vert^2 \f$.  See also
     * GaussianFactorGraph::augmentedJacobian and
     * GaussianFactorGraph::sparseJacobian.
     */
    std::pair<Eigen::MatrixXd,Eigen::VectorXd> jacobian(const Ordering& optionalOrdering) const;
    std::pair<Eigen::MatrixXd,Eigen::VectorXd> jacobian() const;
    /*
      Return a dense \f$ \Lambda \in \mathbb{R}^{n+1 \times n+1} \f$ Hessian
      matrix, augmented with the information vector \f$ \eta \f$.  The
      augmented Hessian is
     \f[ \left[ \begin{array}{ccc}
     \Lambda & \eta \\
     \eta^T & c
     \end{array} \right] \f]
     and the negative log-likelihood is
     \f$ \frac{1}{2} x^T \Lambda x + \eta^T x + c \f$.
     */
    Eigen::MatrixXd augmentedHessian(const Ordering& optionalOrdering) const;
    Eigen::MatrixXd augmentedHessian() const;

    /**
     * Return the dense Hessian \f$ \Lambda \f$ and information vector
     * \f$ \eta \f$, with the noise models baked in. The negative log-likelihood
     * is \frac{1}{2} x^T \Lambda x + \eta^T x + c.  See also
     * GaussianFactorGraph::augmentedHessian.
     */
    std::pair<Eigen::MatrixXd,Eigen::VectorXd> hessian() const;
    std::pair<Eigen::MatrixXd,Eigen::VectorXd> hessian(const Ordering& optionalOrdering) const;

    /** Return only the diagonal of the Hessian A'*A, as a VectorValues */
    virtual std::map<int,Eigen::VectorXd> hessianDiagonal() const;

    /** Return the block diagonal of the Hessian for this factor */
    virtual std::map<int,Eigen::MatrixXd> hessianBlockDiagonal() const;

    /** Solve the factor graph by performing multifrontal variable elimination in COLAMD order using
     *  the dense elimination function specified in \c function (default EliminatePreferCholesky),
     *  followed by back-substitution in the Bayes tree resulting from elimination.  Is equivalent
     *  to calling graph.eliminateMultifrontal()->optimize().*/
    //std::map<int,Eigen::MatrixXd> optimize(const int EliminateKind=1) const;

    std::map<int,Eigen::VectorXd> optimize(Ordering& ordering,
                                           const int EliminateKind=1) const;
    // std::map<int,Eigen::MatrixXd> optimize() const;
    /**
     * Optimize using Eigen's dense Cholesky factorization
     */
    std::map<int,Eigen::VectorXd> optimizeDensely() const;

    /**
     * Compute the gradient of the energy function,
     * \f$ \nabla_{x=x_0} \left\Vert \Sigma^{-1} A x - b \right\Vert^2 \f$,
     * centered around \f$ x = x_0 \f$.
     * The gradient is \f$ A^T(Ax-b) \f$.
     * @param fg The Jacobian factor graph $(A,b)$
     * @param x0 The center about which to compute the gradient
     * @return The gradient as a VectorValues
     */
    std::map<int,Eigen::VectorXd> gradient(const std::map<int,Eigen::VectorXd>& x0) const;

    /**
     * Compute the gradient of the energy function, \f$ \nabla_{x=0} \left\Vert \Sigma^{-1} A x - b
     * \right\Vert^2 \f$, centered around zero. The gradient is \f$ A^T(Ax-b) \f$.
     * @param fg The Jacobian factor graph $(A,b)$
     * @param [output] g A VectorValues to store the gradient, which must be preallocated,
     *        see allocateVectorValues
     * @return The gradient as a VectorValues */
    virtual std::map<int,Eigen::VectorXd> gradientAtZero() const;

    /** Optimize along the gradient direction, with a closed-form computation to perform the line
     *  search.  The gradient is computed about \f$ \delta x=0 \f$.
     *
     *  This function returns \f$ \delta x \f$ that minimizes a reparametrized problem.  The error
     *  function of a GaussianBayesNet is
     *
     *  \f[ f(\delta x) = \frac{1}{2} |R \delta x - d|^2 = \frac{1}{2}d^T d - d^T R \delta x +
     *  \frac{1}{2} \delta x^T R^T R \delta x \f]
     *
     *  with gradient and Hessian
     *
     *  \f[ g(\delta x) = R^T(R\delta x - d), \qquad G(\delta x) = R^T R. \f]
     *
     *  This function performs the line search in the direction of the gradient evaluated at \f$ g =
     *  g(\delta x = 0) \f$ with step size \f$ \alpha \f$ that minimizes \f$ f(\delta x = \alpha g)
     *  \f$:
     *
     *  \f[ f(\alpha) = \frac{1}{2} d^T d + g^T \delta x + \frac{1}{2} \alpha^2 g^T G g \f]
     *
     *  Optimizing by setting the derivative to zero yields \f$ \hat \alpha = (-g^T g) / (g^T G g)
     *  \f$.  For efficiency, this function evaluates the denominator without computing the Hessian
     *  \f$ G \f$, returning
     *
     *  \f[ \delta x = \hat\alpha g = \frac{-g^T g}{(R g)^T(R g)} \f] */
    std::map<int,Eigen::VectorXd> optimizeGradientSearch() const;

    /** x = A'*e */
    std::map<int,Eigen::VectorXd> transposeMultiply(const std::list<Eigen::VectorXd>& e) const;

    /** x += alpha*A'*e */
    void transposeMultiplyAdd(double alpha, const std::list<Eigen::VectorXd>& e, std::map<int,Eigen::VectorXd>& x) const;

    /** return A*x-b */
    std::list<Eigen::VectorXd> gaussianErrors(const std::map<int,Eigen::VectorXd>& x) const;

    ///** return A*x */
    std::list<Eigen::VectorXd> operator*(const std::map<int,Eigen::VectorXd>& x) const;

    ///** y += alpha*A'A*x */
    void multiplyHessianAdd(double alpha, const std::map<int,Eigen::VectorXd>& x,
                            std::map<int,Eigen::VectorXd>& y) const;

    ///** In-place version e <- A*x that overwrites e. */
    void multiplyInPlace(const std::map<int,Eigen::VectorXd>& x, std::list<Eigen::VectorXd>& e) const;

    /** In-place version e <- A*x that takes an iterator. */
    void multiplyInPlace(const std::map<int,Eigen::VectorXd>& x, const std::list<Eigen::VectorXd>::iterator& e) const;

    /// @}


};

/**
 * Evaluates whether linear factors have any constrained noise models
 * @return true if any factor is constrained.
 */
bool hasConstraints(const GaussianFactorGraph& factors);
bool hasConstraints(const std::vector<const RealGaussianFactor*>& factors);

GaussianFactorGraph& getAGFpointer();


std::map<int, std::vector<int>> VariableSlots(const GaussianFactorGraph& factorGraph);
std::map<int, std::vector<int>> VariableSlots(const std::vector<const RealGaussianFactor*>& factorGraph);

std::map<int, std::vector<int>> VariableSlots(std::vector<RealGaussianFactor>& factors1,
 std::vector<RealGaussianFactor>& factors2);

JacobianFactor* convertToJacobianFactorPtr(RealGaussianFactor* gf);
const JacobianFactor* convertToJacobianFactorPtr(const RealGaussianFactor* gf);

const JacobianFactor* convertToJacobianFactorConstPtr(const RealGaussianFactor* gf);
#endif // GAUSSIANFACTORPOINTERGRAPH_H_INCLUDED
