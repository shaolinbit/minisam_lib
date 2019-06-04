#ifndef BAYESTREEPOINTER_H_INCLUDED
#define BAYESTREEPOINTER_H_INCLUDED

/**
 * @file    BayesTree.h
 * @brief   Bayes Tree is a tree of cliques of a Bayes Chain
 * @author
 */

#include "../linear/GaussianFactorGraph.h"
#include "../inference/BayesTreeCliqueBasePointer.h"
#include "../linear/GaussianConditional.h"
#include "../linear/GaussianBayesNetPointer.h"

#include <string>


class GaussianFactorGraph;
/* ************************************************************************* */
/** clique statistics */
struct BayesTreeCliquePointerStats
{
    double avgConditionalSize;
    int maxConditionalSize;
    double avgSeparatorSize;
    int maxSeparatorSize;
    //void print(const std::string& s = "") const ;
};

/** store all the sizes  */
struct BayesTreeCliquePointerData
{
    std::vector<int> conditionalSizes;
    std::vector<int> separatorSizes;
    BayesTreeCliquePointerStats getStats() const;
};

/* ************************************************************************* */
/**
 * Bayes tree
 * @tparam CONDITIONAL The type of the conditional densities, i.e. the type of node in the underlying Bayes chain,
 * which could be a ConditionalProbabilityTable, a GaussianConditional, or a SymbolicConditional.
 * @tparam CLIQUE The type of the clique data structure, defaults to BayesTreeClique, normally do not change this
 * as it is only used when developing special versions of BayesTree, e.g. for ISAM2.
 *
 * \addtogroup Multifrontal
 * \nosubgrouping
 */
class BayesTreePointer
{
public:
    /** Map from indices to Clique */
    std::map<int,BayesTreeCliqueBasePointer*> nodesbtc;
    std::vector<BayesTreeCliqueBasePointer*> roots_;

    /// @name Standard Constructors
    /// @{
public:
    /** Create an empty Bayes Tree */
    BayesTreePointer();

    virtual ~BayesTreePointer();

    /** Copy constructor */
    BayesTreePointer(BayesTreePointer& other);

    // BayesTree(const BayesTree& other);

    BayesTreePointer(const std::map<int,BayesTreeCliqueBasePointer*>& nbtc,const std::vector<BayesTreeCliqueBasePointer*>& nroots):
        nodesbtc(nbtc),roots_(nroots) {}

    BayesTreePointer(const BayesTreePointer& other) = default;

    /** Assignment operator */

    BayesTreePointer operator=(const BayesTreePointer& other);

    /// @name Testable
    /// @{

public:

    /// @name Standard Interface
    /// @{

    /** number of cliques */
    int size() const;

    /** Check if there are any cliques in the tree */
    inline bool empty() const
    {
        return (nodesbtc.empty()&&roots_.empty());
    }


    /** Access node by variable */
    BayesTreeCliqueBasePointer* operator[](int j) const;
    BayesTreeCliqueBasePointer* operator[](int j);
    /** return root cliques */
    const std::vector<BayesTreeCliqueBasePointer* > roots() const
    {
        return roots_;
    }

    /** alternate syntax for matlab: find the clique that contains the variable with Key j */
    BayesTreeCliqueBasePointer* clique(int j);

    /** Gather data on all cliques */
    BayesTreeCliquePointerData getCliqueData() const;

    /** Collect number of cliques with cached separator marginals */
      int numCachedSeparatorMarginals() const;

    /** Return marginal on any variable.  Note that this actually returns a conditional, for which a
     *  solution may be directly obtained by calling .solve() on the returned object.
     *  Alternatively, it may be directly used as its factor base class.  For example, for Gaussian
     *  systems, this returns a GaussianConditional, which inherits from JacobianFactor and
     *  GaussianFactor. */
    GaussianConditional* marginalFactor(int j,int Eliminatefunction);

    /**
     * return joint on two variables
     * Limitation: can only calculate joint if cliques are disjoint or one of them is root
     */
    GaussianFactorGraph* joint(int j1, int j2, const int Eliminatefunction);

    /**marginalFactor
     * return joint on two variables as a BayesNet
     * Limitation: can only calculate joint if cliques are disjoint or one of them is root
     */
    GaussianBayesNetPointer* jointBayesNet(int j1, int j2, const int Eliminatefunction);

    /**
     * Read only with side effects
     */

    /** saves the Tree to a text file in GraphViz format */
    void saveGraph(std::ostream& stm) const;



    /// @}
    /// @name Advanced Interface
    /// @{

    /**
     * Find parent clique of a conditional.  It will look at all parents and
     * return the one with the lowest index in the ordering.
     */
    int findParentClique(const GaussianConditional& parents) const;

    /** Remove all nodes */
    void clear();

    /** Clear all shortcut caches - use before timing on marginal calculation to avoid residual cache data */
    void deleteCachedShortcuts();

    /**
     * Remove path from clique to root and return that path as factors
     * plus a list of orphaned subtree roots. Used in removeTop below.
     */
    void removePath(BayesTreeCliqueBasePointer* clique, GaussianBayesNetPointer* bn,
                                std::list<BayesTreeCliqueBasePointer*> *orphans);
    void fakeremovePath(BayesTreeCliqueBasePointer* clique);

    /**
     * Given a list of indices, turn "contaminated" part of the tree back into a factor graph.
     * Factors and orphans are added to the in/out arguments.
     */
    void removeTop(const std::vector<int>& keys, GaussianBayesNetPointer* bn, std::list<BayesTreeCliqueBasePointer*> *orphans);

    /**
     * Remove the requested subtree. */
    std::list<BayesTreeCliqueBasePointer*> removeSubtree(BayesTreeCliqueBasePointer* subtree);

    /** Insert a new subtree with known parent clique.  This function does not check that the
     *  specified parent is the correct parent.  This function updates all of the internal data
     *  structures associated with adding a subtree, such as populating the nodes index. */
    void insertRoot(BayesTreeCliqueBasePointer* subtree);

    /** add a clique (top down) */
  //  void addClique(BayesTreeCliqueBase& clique, BayesTreeCliqueBase& parent_clique);

    /** Add all cliques in this BayesTree to the specified factor graph */
    void addFactorsToGraph(GaussianFactorGraph& graph) const;

    void ConvertToGraph(GaussianFactorGraph& nb) const;

    void addconditionaltograph(GaussianFactorGraph& nb,const BayesTreeCliqueBasePointer* cc) const;

    /** Recursively optimize the BayesTree to produce a vector solution. */
    std::map<int,Eigen::VectorXd> optimize() const;


    /**
         * Optimize along the gradient direction, with a closed-form computation to perform the line
         * search.  The gradient is computed about \f$ \delta x=0 \f$.
         *
         * This function returns \f$ \delta x \f$ that minimizes a reparametrized problem.  The error
         * function of a GaussianBayesNet is
         *
         * \f[ f(\delta x) = \frac{1}{2} |R \delta x - d|^2 = \frac{1}{2}d^T d - d^T R \delta x +
         * \frac{1}{2} \delta x^T R^T R \delta x \f]
         *
         * with gradient and Hessian
         *
         * \f[ g(\delta x) = R^T(R\delta x - d), \qquad G(\delta x) = R^T R. \f]
         *
         * This function performs the line search in the direction of the gradient evaluated at \f$ g =
         * g(\delta x = 0) \f$ with step size \f$ \alpha \f$ that minimizes \f$ f(\delta x = \alpha g)
         * \f$:
         *
         * \f[ f(\alpha) = \frac{1}{2} d^T d + g^T \delta x + \frac{1}{2} \alpha^2 g^T G g \f]
         *
         * Optimizing by setting the derivative to zero yields \f$ \hat \alpha = (-g^T g) / (g^T G g)
         * \f$.  For efficiency, this function evaluates the denominator without computing the Hessian
         * \f$ G \f$, returning
         *
         * \f[ \delta x = \hat\alpha g = \frac{-g^T g}{(R g)^T(R g)} \f] */
    std::map<int,Eigen::VectorXd> optimizeGradientSearch();

    /** Compute the gradient of the energy function, \f$ \nabla_{x=x_0} \left\Vert \Sigma^{-1} R x -
     * d \right\Vert^2 \f$, centered around \f$ x = x_0 \f$. The gradient is \f$ R^T(Rx-d) \f$.
     *
     * @param x0 The center about which to compute the gradient
     * @return The gradient as a VectorValues */
    std::map<int,Eigen::VectorXd> gradient(const std::map<int,Eigen::VectorXd>& x0);

    /** Compute the gradient of the energy function, \f$ \nabla_{x=0} \left\Vert \Sigma^{-1} R x - d
     * \right\Vert^2 \f$, centered around zero. The gradient about zero is \f$ -R^T d \f$.  See also
     * gradient(const GaussianBayesNet&, const VectorValues&).
     *
     * @return A VectorValues storing the gradient. */
    std::map<int,Eigen::VectorXd> gradientAtZero();

    /** Mahalanobis norm error. */
    double error(const std::map<int,Eigen::VectorXd>& x) const;

    /** Computes the determinant of a GassianBayesTree, as if the Bayes tree is reorganized into a
     * matrix. A GassianBayesTree is equivalent to an upper triangular matrix, and for an upper
     * triangular matrix determinant is the product of the diagonal elements. Instead of actually
     * multiplying we add the logarithms of the diagonal elements and take the exponent at the end
     * because this is more numerically stable. */
    double determinant() const;

    /** Computes the determinant of a GassianBayesTree, as if the Bayes tree is reorganized into a
     * matrix. A GassianBayesTree is equivalent to an upper triangular matrix, and for an upper
     * triangular matrix determinant is the product of the diagonal elements. Instead of actually
     * multiplying we add the logarithms of the diagonal elements and take the exponent at the end
     * because this is more numerically stable. */
    double logDeterminant() const;

    /** Return the marginal on the requested variable as a covariance matrix.  See also
    *   marginalFactor(). */
    Eigen::MatrixXd marginalCovariance(int key);




protected:

    /** private helper method for saving the Tree to a text file in GraphViz format */
    void saveGraph(std::ostream &s,const BayesTreeCliqueBasePointer* clique,
                   int parentnum = 0) const;

    /** Gather data on a single clique */
    void getCliqueData(BayesTreeCliquePointerData& data,const BayesTreeCliqueBasePointer* clique) const;

    /** remove a clique: warning, can result in a forest */
    std::vector<int> removeClique(BayesTreeCliqueBasePointer* clique);

    /** Fill the nodes index for a subtree */
     void fillNodesIndex(BayesTreeCliqueBasePointer* subtree);

    // Friend JunctionTree because it directly fills roots and nodes index.
    friend class EliminatableClusterTree;

}; // BayesTree

/* ************************************************************************* */
///not needed I strongly believed this class is not necessary.

//template<class CLIQUE>
class BayesTreeOrphanWrapperPointer : public JacobianFactor
{
public:
    BayesTreeCliqueBasePointer* clique_;
public:
    BayesTreeOrphanWrapperPointer(BayesTreeCliqueBasePointer* clique,const std::vector<int>& inkeys):JacobianFactor(JacobianFactor())
    {
        // Store parent keys in our base type factor so that eliminating those parent keys will pull
        // this subtree into the elimination.
        this->keys_=inkeys;
        clique_=clique;
        this->iswrapper_=true;
    }
};
struct BayesTreeTraversalNodeINTPointer
{
    bool expanded;
    BayesTreeCliqueBasePointer* treeNode;
    int parentData;
    //typename FastList<DATA>::iterator dataPointer;
    //std::list<int>::iterator
    BayesTreeTraversalNodeINTPointer(BayesTreeCliqueBasePointer* _treeNode, int _parentData) :
        expanded(false), treeNode(_treeNode), parentData(_parentData)
    {
    }
};
struct BayesTreeTraversalNodeDPointer
{
    bool expanded;
    BayesTreeCliqueBasePointer* treeNode;
    double parentData;
    //typename FastList<DATA>::iterator dataPointer;
    //std::list<int>::iterator
    BayesTreeTraversalNodeDPointer(BayesTreeCliqueBasePointer* _treeNode, double _parentData) :
        expanded(false), treeNode(_treeNode), parentData(_parentData)
    {
    }
};


struct OptimizeDataPointer
{
     OptimizeDataPointer* parentData_;
    std::map<int, std::map<int,Eigen::VectorXd>::const_iterator> cliqueResults;
    OptimizeDataPointer(OptimizeDataPointer* bparentData,std::map<int, std::map<int,Eigen::VectorXd>::const_iterator>
                        bcliqueResults):parentData_(bparentData),cliqueResults(bcliqueResults) {}
};

struct BT_LATraversalNodePointer
{
    bool expanded;
    BayesTreeCliqueBasePointer* treeNode;
    OptimizeDataPointer* parentData;
    //typename FastList<DATA>::iterator dataPointer;
    std::list<OptimizeDataPointer*>::iterator dataPointer;
    BT_LATraversalNodePointer(BayesTreeCliqueBasePointer* _treeNode, OptimizeDataPointer* _parentData) :
        expanded(false), treeNode(_treeNode), parentData(_parentData)
    {
    }
};

OptimizeDataPointer* OptimizeCliquePointer(std::map<int,Eigen::VectorXd>& collectedResult,BayesTreeCliqueBasePointer* clique,
        OptimizeDataPointer* parentData);

void BayesTreeDepthFirstForestINTPointer(const BayesTreePointer& forest, int* rootData,GaussianFactorGraph& graph);

void GBayesTreeDepthFirstForestdoublePointer(const BayesTreePointer& forest, double* rootData);


//void BayesTreeDepthFirstForestBTCB(BayesTree& forest, BayesTreeCliqueBase* rootData);


/** Compute the marginal of the requested variables and return the result as a Bayes net.
    *  @param variables Determines the variables whose marginal to compute, if provided as an
    *         Ordering they will be ordered in the returned BayesNet as specified, and if provided
    *         as a vector<Key> they will be ordered using constrained COLAMD.
    *  @param marginalizedVariableOrdering Optional ordering for the variables being marginalized
    *         out, i.e. all variables not in \c variables.  If this is boost::none, the ordering
    *         will be computed with COLAMD.
    *  @param function Optional dense elimination function, if not provided the default will be
    *         used.
    *  @param variableIndex Optional pre-computed VariableIndex for the factor graph, if not
    *         provided one will be computed.
   GaussianBayesNet marginalMultifrontalBayesNet(const GaussianFactorGraph& gf,
     Ordering& variables,
     const int Eliminatefunction);*/

// BayesTreeCliqueBase
//BayesTreeCloneForestVisitorPre(const BayesTreeCliqueBase& node, BayesTreeCliqueBase* parentPointer);

BayesTreePointer* makeNullBayesTreePointer();

std::map<int,Eigen::VectorXd> optimizeBayesTreePointer(const BayesTreePointer& bayesTree);

double internallogDeterminantPointer(BayesTreeCliqueBasePointer* clique, double* parentSum);



#endif // BAYESTREEPOINTER_H_INCLUDED
