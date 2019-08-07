#ifndef BAYESTREE_H_INCLUDED
#define BAYESTREE_H_INCLUDED

/**
 * @file    BayesTree.h
 * @brief   Bayes Tree is a tree of cliques of a Bayes Chain
 * @author  Frank Dellaert
 */

#include "../linear/GaussianFactorGraph.h"
#include "../inference/BayesTreeCliqueBase.h"
#include "../linear/GaussianConditional.h"
#include "../linear/GaussianBayesNet.h"

#include <string>

namespace minisam
{
class GaussianFactorGraph;

/* ************************************************************************* */
/** clique statistics */
struct BayesTreeCliqueStats
{
    double avgConditionalSize;
    int maxConditionalSize;
    double avgSeparatorSize;
    int maxSeparatorSize;
};

/** store all the sizes  */
struct BayesTreeCliqueData
{
    std::vector<int> conditionalSizes;
    std::vector<int> separatorSizes;
    BayesTreeCliqueStats getStats() const;
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
class BayesTree
{
public:
    /** Map from indices to Clique */
    std::map<int,BayesTreeCliqueBase*>* nodesbtc;
    std::vector<BayesTreeCliqueBase*>* roots_;

    /// @name Standard Constructors
    /// @{
public:
    /** Create an empty Bayes Tree */
    BayesTree();

    virtual ~BayesTree();

    /** Copy constructor */
    BayesTree(BayesTree& other);

    BayesTree(const BayesTree& other) = default;

    /** Assignment operator */

    BayesTree operator=(const BayesTree& other);

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
        return ((nodesbtc->empty())&&(roots_->empty()));
    }


    /** Access node by variable */
    BayesTreeCliqueBase* operator[](int j) const;
    BayesTreeCliqueBase* operator[](int j);
    /** return root cliques */
    std::vector<BayesTreeCliqueBase* >* roots() const;

    /** alternate syntax for matlab: find the clique that contains the variable with Key j */
    BayesTreeCliqueBase* clique(int j);

    /** Gather data on all cliques */
    BayesTreeCliqueData getCliqueData() const;

    /** Collect number of cliques with cached separator marginals */
    int numCachedSeparatorMarginals() const;

    /** Return marginal on any variable.  Note that this actually returns a conditional, for which a
     *  solution may be directly obtained by calling .solve() on the returned object.
     *  Alternatively, it may be directly used as its factor base class.  For example, for Gaussian
     *  systems, this returns a GaussianConditional, which inherits from JacobianFactor and
     *  GaussianFactor. */
    std::pair<GaussianConditional*,GaussianBayesNet*> marginalFactor(int j,int Eliminatefunction);

    /**
     * return joint on two variables
     * Limitation: can only calculate joint if cliques are disjoint or one of them is root
     */
    GaussianFactorGraph* joint(int j1, int j2, const int Eliminatefunction);

    /**marginalFactor
     * return joint on two variables as a BayesNet
     * Limitation: can only calculate joint if cliques are disjoint or one of them is root
     */
    GaussianBayesNet* jointBayesNet(int j1, int j2, const int Eliminatefunction);

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
    void clearallcliques();
    void clearall();

    /** Clear all shortcut caches - use before timing on marginal calculation to avoid residual cache data */
    void deleteCachedShortcuts();

    /**
     * Remove path from clique to root and return that path as factors
     * plus a list of orphaned subtree roots. Used in removeTop below.
     */
    void removePath(BayesTreeCliqueBase* clique, FactorGraph<GaussianConditional>* bn,
                    std::list<BayesTreeCliqueBase*> *orphans);
    /**
     * Given a list of indices, turn "contaminated" part of the tree back into a factor graph.
     * Factors and orphans are added to the in/out arguments.
     */
    void removeTop(const std::vector<int>& keys, GaussianBayesNet* bn, std::list<BayesTreeCliqueBase*> *orphans);


    /** Add all cliques in this BayesTree to the specified factor graph */
    void addFactorsToGraph(GaussianFactorGraph& graph) const;

    void ConvertToGraph(GaussianFactorGraph& nb) const;

    void addconditionaltograph(GaussianFactorGraph& nb,const BayesTreeCliqueBase* cc) const;

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




protected:

    /** private helper method for saving the Tree to a text file in GraphViz format */
    void saveGraph(std::ostream &s,const BayesTreeCliqueBase* clique,
                   int parentnum = 0) const;

    /** Gather data on a single clique */
    void getCliqueData(BayesTreeCliqueData& data,const BayesTreeCliqueBase* clique) const;

    /** remove a clique: warning, can result in a forest */
    int removeClique(BayesTreeCliqueBase* clique);

    /** Fill the nodes index for a subtree */
    void fillNodesIndex(BayesTreeCliqueBase* subtree);

    // Friend JunctionTree because it directly fills roots and nodes index.
    friend class EliminatableClusterTree;

}; // BayesTree

/* ************************************************************************* */

class BayesTreeOrphanWrapper : public JacobianFactor
{
public:
    BayesTreeCliqueBase* clique_;
public:
    BayesTreeOrphanWrapper(BayesTreeCliqueBase* clique,const std::vector<int>& inkeys):JacobianFactor(JacobianFactor())
    {
        // Store parent keys in our base type factor so that eliminating those parent keys will pull
        // this subtree into the elimination.
        this->keys_=inkeys;
        clique_=clique;
        this->iswrapper_=true;
    }
};
struct BayesTreeTraversalNodeINT
{
    bool expanded;
    BayesTreeCliqueBase* treeNode;
    int parentData;
    BayesTreeTraversalNodeINT(BayesTreeCliqueBase* _treeNode, int _parentData) :
        expanded(false), treeNode(_treeNode), parentData(_parentData)
    {
    }
};
struct BayesTreeTraversalNodeDPointer
{
    bool expanded;
    BayesTreeCliqueBase* treeNode;
    double parentData;
    BayesTreeTraversalNodeDPointer(BayesTreeCliqueBase* _treeNode, double _parentData) :
        expanded(false), treeNode(_treeNode), parentData(_parentData)
    {
    }
};


struct OptimizeData
{
    OptimizeData* parentData_;
    std::map<int, std::map<int,Eigen::VectorXd>::const_iterator> cliqueResults;
    OptimizeData(OptimizeData* bparentData,std::map<int, std::map<int,Eigen::VectorXd>::const_iterator>
                 bcliqueResults):parentData_(bparentData),cliqueResults(bcliqueResults) {}
    OptimizeData(const OptimizeData& opti):parentData_(opti.parentData_),cliqueResults(opti.cliqueResults) {}

};

struct BT_LATraversalNode
{
    bool expanded;
    const BayesTreeCliqueBase* treeNode;
    OptimizeData* parentData;
    std::list<OptimizeData>::iterator dataPointer;
    BT_LATraversalNode(const BayesTreeCliqueBase* _treeNode, OptimizeData* _parentData) :
        expanded(false), treeNode(_treeNode), parentData(_parentData)
    {
    }
};

OptimizeData OptimizeClique(std::map<int,Eigen::VectorXd>& collectedResult,const BayesTreeCliqueBase* clique,
                            OptimizeData* parentData);

void BayesTreeDepthFirstForestINT(const BayesTree& forest, int* rootData,GaussianFactorGraph& graph);

void GBayesTreeDepthFirstForestdouble(const BayesTree& forest, double* rootData);


std::map<int,Eigen::VectorXd> optimizeBayesTree(const BayesTree& bayesTree);

double internallogDeterminant(BayesTreeCliqueBase* clique, double* parentSum);
};

#endif // BAYESTREE_H_INCLUDED
