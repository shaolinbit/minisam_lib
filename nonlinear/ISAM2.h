#ifndef ISAM2POINTER_H_INCLUDED
#define ISAM2POINTER_H_INCLUDED


#include "../nonlinear/ISAM2-implPointer.h"
#include "../nonlinear/NonlinearFactorGraph.h"
#include "../nonlinear/DoglegOptimizerImpl.h"
//#include "../inference/BayesTreePointer.h"
#include "../nonlinear/ISAM2Params.h"
#include "../inference/VariableIndex.h"
#include "../nonlinear/ISAM2CliquePointer.h"

/**
 * @file    ISAM2.h
 * @brief   Incremental update functionality (ISAM2) for BayesTree, with fluid relinearization.
 * @author
 */


/**
 * @addtogroup ISAM2
 * Implementation of the full ISAM2 algorithm for incremental nonlinear optimization.
 *
 * The typical cycle of using this class to create an instance by providing ISAM2Params
 * to the constructor, then add measurements and variables as they arrive using the update()
 * method.  At any time, calculateEstimate() may be called to obtain the current
 * estimate of all variables.
 *
 */
class ISAM2Data
{
public:
     /** The current linearization point */
    std::map<int,Eigen::VectorXd> theta_;

#ifdef GMF_Using_Pose3
std::map<int,Pose3> thetaPose_;
#else
    std::map<int,Pose2> thetaPose_;
#endif // GMF_Using_Pose3

  mutable std::map<int,Eigen::VectorXd> delta_;

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
      /// Access the current linearization point
    const std::map<int,Eigen::VectorXd>& getLinearizationPoint() const
    {
        return theta_;
    }

#ifdef GMF_Using_Pose3
    const std::map<int,Pose3>& getLinearizationPose() const
    {
        return thetaPose_;
    }
#else
    const std::map<int,Pose2>& getLinearizationPose() const
    {
        return thetaPose_;
    }
#endif // GMF_Using_Pose3

void clearfactors();
void clearpose();
 /** Access the set of nonlinear factors*/
    const NonlinearFactorGraph getFactorsUnsafe() const
    {
        return nonlinearFactors_;
    }

    /** Access the nonlinear variable index */
    const VariableIndex getVariableIndex() const
    {
        return variableIndex_;
    }

    /** Access the nonlinear variable index */
    const std::set<int> getFixedVariables() const
    {
        return fixedVariables_;
    }

};


class  ISAM2
{
public:
    /** Map from indices to Clique */
   // std::map<int,int> nodeskey_index;
    /**CliqueBase vector**/
    //std::vector<ISAM2CliquePointer*> nodesbtc;
    /** Root cliques */
   // std::vector<int> roots_;
   // int rootcliqueindex;
   std::map<int,ISAM2CliquePointer*>* nodesbtc;
   std::vector<ISAM2CliquePointer*>* roots_;
   //std::vector<ISAM2CliquePointer*> nodesvector;


//protected:

    /** The current linearization point */
   // std::map<int,Eigen::VectorXd> theta_;

////#ifdef GMF_Using_Pose3
  //  std::map<int,Pose3*> thetaPose_;
//#else
//    std::map<int,Pose2*> thetaPose_;
//#endif // GMF_Using_Pose3
    /** The linear delta from the last linear solution, an update to the estimate in theta
      *
      * This is \c mutable because it is a "cached" variable - it is not updated
      * until either requested with getDelta() or calculateEstimate(), or needed
      * during update() to evaluate whether to relinearize variables.
      */
   // mutable std::map<int,Eigen::VectorXd> delta_;

   // mutable std::map<int,Eigen::VectorXd> deltaNewton_; // Only used when using Dogleg - stores the Gauss-Newton update
  //  mutable std::map<int,Eigen::VectorXd> RgProd_; // Only used when using Dogleg - stores R*g and is updated incrementally

    /** The current parameters */
    ISAM2Params params_;
   // NonlinearFactorGraph nonlinearFactors_;
    /** The current linear factors, which are only updated as needed */
   // mutable GaussianFactorGraph linearFactors_;
protected:
    /** VariableIndex lets us look up factors by involved variable and keeps track of dimensions */
   // VariableIndex variableIndex_;
    /** All original nonlinear factors are stored here to use during relinearization */


    /** A cumulative mask for the variables that were replaced and have not yet
     * been updated in the linear solution delta_, this is only used internally,
     * delta will always be updated if necessary when requested with getDelta()
     * or calculateEstimate().
     *
     * This is \c mutable because it is used internally to not update delta_
     * until it is needed.
     */
 //   mutable std::set<int> deltaReplacedMask_; // TODO: Make sure accessed in the right way





    /** The current Dogleg Delta (trust region radius) */
    mutable double doglegDelta_;

    /** Set of variables that are involved with linear factors from marginalized
     * variables and thus cannot have their linearization points changed. */
    //std::set<int> fixedVariables_;

    int update_count_; ///< Counter incremented every update(), used to determine periodic relinearization

public:

    /** Create an empty ISAM2 instance */
    ISAM2(const ISAM2Params& params);

    /** Create an empty ISAM2 instance using the default set of parameters (see ISAM2Params) */
    ISAM2();

    /** default virtual destructor */
    ~ISAM2();


    ISAM2(const ISAM2& other);

    ISAM2& operator=(const ISAM2& other);

    // ISAM2(const ISAM2&) = default;

    int updatetimefortune;

    /**
     * Add new factors, updating the solution and relinearizing as needed.
     *
     * Optionally, this function remove existing factors from the system to enable
     * behaviors such as swapping existing factors with new ones.
     *
     * Add new measurements, and optionally new variables, to the current system.
     * This runs a full step of the ISAM2 algorithm, relinearizing and updating
     * the solution as needed, according to the wildfire and relinearize
     * thresholds.
     *
     * @param newFactors The new factors to be added to the system
     * @param newTheta Initialization points for new variables to be added to the system.
     * You must include here all new variables occuring in newFactors (which were not already
     * in the system).  There must not be any variables here that do not occur in newFactors,
     * and additionally, variables that were already in the system must not be included here.
     * @param removeFactorIndices Indices of factors to remove from system
     * @param force_relinearize Relinearize any variables whose delta magnitude is sufficiently
     * large (Params::relinearizeThreshold), regardless of the relinearization interval
     * (Params::relinearizeSkip).
     * @param constrainedKeys is an optional map of keys to group labels, such that a variable can
     * be constrained to a particular grouping in the BayesTree
     * @param noRelinKeys is an optional set of nonlinear keys that iSAM2 will hold at a constant linearization
     * point, regardless of the size of the linear delta
     * @param extraReelimKeys is an optional set of nonlinear keys that iSAM2 will re-eliminate, regardless
     * of the size of the linear delta. This allows the provided keys to be reordered.
     * @return An ISAM2Result struct containing information about the update

    ISAM2Result update(NonlinearFactorGraph& newFactors,
        const std::map<int,Eigen::VectorXd>& newTheta,
        std::vector<int>* removeFactorIndices=NULL,
         std::map<int,int> *constrainedKeys=NULL,//const boost::optional<FastMap<Key,int> >& constrainedKeys = boost::none,
        std::list<int> *noRelinKeys=NULL,//const boost::optional<FastList<Key> >& noRelinKeys = boost::none,
        std::list<int> *extraReelimKeys=NULL,//const boost::optional<FastList<Key> >& extraReelimKeys = boost::none,
        bool force_relinearize = false);*/
#ifdef GMF_Using_Pose3
    ISAM2Result update(NonlinearFactorGraph& newFactors,
                       const std::map<int,Eigen::VectorXd>& newTheta,
                       const std::map<int,Pose3>& newposetheta,
                       ISAM2Data& isam2data,
                       std::vector<int>* removeFactorIndices=NULL,
                       std::map<int,int> *constrainedKeys=NULL,//const boost::optional<FastMap<Key,int> >& constrainedKeys = boost::none,
                       std::list<int> *noRelinKeys=NULL,//const boost::optional<FastList<Key> >& noRelinKeys = boost::none,
                       std::list<int> *extraReelimKeys=NULL,//const boost::optional<FastList<Key> >& extraReelimKeys = boost::none,
                       bool force_relinearize = false);
#else
    ISAM2Result update(NonlinearFactorGraph& newFactors,
                       const std::map<int,Eigen::VectorXd>& newTheta,
                       const std::map<int,Pose2>& newposetheta,
                       ISAM2Data& isam2data,
                       std::vector<int>* removeFactorIndices=NULL,
                       std::map<int,int> *constrainedKeys=NULL,//const boost::optional<FastMap<Key,int> >& constrainedKeys = boost::none,
                       std::list<int> *noRelinKeys=NULL,//const boost::optional<FastList<Key> >& noRelinKeys = boost::none,
                       std::list<int> *extraReelimKeys=NULL,//const boost::optional<FastList<Key> >& extraReelimKeys = boost::none,
                       bool force_relinearize = false);
#endif // GMF_Using_Pose3
    /** Marginalize out variables listed in leafKeys.  These keys must be leaves
     * in the BayesTree.  Throws MarginalizeNonleafException if non-leaves are
     * requested to be marginalized.  Marginalization leaves a linear
     * approximation of the marginal in the system, and the linearization points
     * of any variables involved in this linear marginal become fixed.  The set
     * fixed variables will include any key involved with the marginalized variables
     * in the original factors, and possibly additional ones due to fill-in.
     *
     * If provided, 'marginalFactorsIndices' will be augmented with the factor graph
     * indices of the marginal factors added during the 'marginalizeLeaves' call
     *
     * If provided, 'deletedFactorsIndices' will be augmented with the factor graph
     * indices of any factor that was removed during the 'marginalizeLeaves' call
     */
    void marginalizeLeaves(const std::list<int>& leafKeys,ISAM2Data& isam2data,
                           std::vector<int> *marginalFactorsIndices=NULL, //boost::optional<FactorIndices&> marginalFactorsIndices = boost::none,
                           std::vector<int> *deletedFactorsIndices=NULL);//boost::optional<FactorIndices&> deletedFactorsIndices = boost::none);
/*
    /// Access the current linearization point
    const std::map<int,Eigen::VectorXd>& getLinearizationPoint() const
    {
        return theta_;
    }

#ifdef GMF_Using_Pose3
    const std::map<int,Pose3*>& getLinearizationPose() const
    {
        return thetaPose_;
    }
#else
    const std::map<int,Pose2*>& getLinearizationPose() const
    {
        return thetaPose_;
    }
#endif // GMF_Using_Pose3
*/


    /** Compute an estimate from the incomplete linear delta computed during the last update.
     * This delta is incomplete because it was not updated below wildfire_threshold.  If only
     * a single variable is needed, it is faster to call calculateEstimate(const KEY&).
     */
    std::map<int,Eigen::VectorXd> calculateEstimate(ISAM2Data& isam2data);

#ifdef GMF_Using_Pose3
    std::map<int,Eigen::VectorXd> calculateEstimate(ISAM2Data& isam2data,std::map<int,Pose3>* pose3lin);
#else
    std::map<int,Eigen::VectorXd> calculateEstimate(ISAM2Data& isam2data,std::map<int,Pose2>* pose2lin);
#endif

    /** Compute an estimate for a single variable using its incomplete linear delta computed
     * during the last update.  This is faster than calling the no-argument version of
     * calculateEstimate, which operates on all variables.
     * @param key
     * @return
     */
    //template<class VALUE>
    Eigen::VectorXd calculateEstimate(ISAM2Data& isam2data,int key);
#ifdef GMF_Using_Pose3
    Pose3 calculateEstimatePose(ISAM2Data& isam2data,int key);
#else
    Pose2 calculateEstimatePose(ISAM2Data& isam2data,int key);
#endif // GMF_Using_Pose3

    /** Compute an estimate for a single variable using its incomplete linear delta computed
     * during the last update.  This is faster than calling the no-argument version of
     * calculateEstimate, which operates on all variables.  This is a non-templated version
     * that returns a Value base class for use with the MATLAB wrapper.
     * @param key
     * @return
     */

    /** Return marginal on any variable as a covariance matrix */
    Eigen::MatrixXd marginalCovariance(int key);


    /** Compute an estimate using a complete delta computed by a full back-substitution.
     */
    std::map<int,Eigen::VectorXd> calculateBestEstimate() const;

    /** Access the current delta, computed during the last call to update */
    std::map<int,Eigen::VectorXd>& getDelta(ISAM2Data& isam2data);

    /** Compute the linear error */
    double error(const std::map<int,Eigen::VectorXd>& x) const;

    /** Access the set of nonlinear factors
    const NonlinearFactorGraph getFactorsUnsafe() const
    {
        return nonlinearFactors_;
    }

    /** Access the nonlinear variable index
    const VariableIndex getVariableIndex() const
    {
        return variableIndex_;
    }

    /** Access the nonlinear variable index
    const std::set<int> getFixedVariables() const
    {
        return fixedVariables_;
    } */

    int lastAffectedVariableCount;
    int lastAffectedFactorCount;
    int lastAffectedCliqueCount;
    int lastAffectedMarkedCount;
    mutable int lastBacksubVariableCount;
    int lastNnzTop;

    const ISAM2Params params() const
    {
        return params_;
    }

    /** Compute the gradient of the energy function, \f$ \nabla_{x=0} \left\Vert \Sigma^{-1} R x - d
     * \right\Vert^2 \f$, centered around zero. The gradient about zero is \f$ -R^T d \f$.  See also
     * gradient(const GaussianBayesNet&, const VectorValues&).
     *
     * @return A VectorValues storing the gradient.
     */
    //std::map<int,Eigen::VectorXd> gradientAtZero() const;



    /// @}


    //functions from Bayestree;

    /** number of cliques */
    int size() const;

    /** Check if there are any cliques in the tree */
    inline bool empty() const
    {
        return (roots()->empty()&&nodesbtc->empty());
    }

    /** return nodes
    const std::vector<int>  getnodeskey() const
    {
        return nodeskey;
    }*/

    /** Access node by variable */
    ISAM2CliquePointer* operator[](int j) const;
    ISAM2CliquePointer* operator[](int j);

    /** return root cliques */
    std::vector<ISAM2CliquePointer*>* roots() const;

    /** alternate syntax for matlab: find the clique that contains the variable with Key j */
    ISAM2CliquePointer* clique(int j);

    /** Gather data on all cliques */
    //  BayesTreeCliqueData getCliqueData() const;

    /** Collect number of cliques with cached separator marginals */
//    int numCachedSeparatorMarginals() const;

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
    void  ConvertToGraph(GaussianFactorGraph& nb);
    void addconditionaltograph(GaussianFactorGraph& nb,const ISAM2CliquePointer* cc);

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
    // template<class CONTAINER>
    int findParentClique(const GaussianConditional& parents) const;

    /** Remove all nodes */
    void clear();
    //void clearfactors();
    void clearcliques();
    void clearallcliques();
    void clearall();

    /** Clear all shortcut caches - use before timing on marginal calculation to avoid residual cache data */
    void deleteCachedShortcuts();

    /**
     * Remove path from clique to root and return that path as factors
     * plus a list of orphaned subtree roots. Used in removeTop below.
     */
    void removePath(ISAM2CliquePointer* clique, GaussianBayesNetPointer* bn,
                    std::list<ISAM2CliquePointer*> *orphans);

    /**
     * Given a list of indices, turn "contaminated" part of the tree back into a factor graph.
     * Factors and orphans are added to the in/out arguments.
     */
    void removeTop(const std::vector<int>& keys, GaussianBayesNetPointer* bn, std::list<ISAM2CliquePointer*> *orphans);

    /**
     * Remove the requested subtree. */
    std::list<ISAM2CliquePointer*> removeSubtree(ISAM2CliquePointer* subtree);

    /** Insert a new subtree with known parent clique.  This function does not check that the
     *  specified parent is the correct parent.  This function updates all of the internal data
     *  structures associated with adding a subtree, such as populating the nodes index. */
    // void insertRoot(const ISAM2Clique& subtree);

    /** add a clique (top down) */
    // void addClique(ISAM2Clique& clique, ISAM2Clique& parent_clique);

    /** Add all cliques in this BayesTree to the specified factor graph */
    //  void addFactorsToGraph(GaussianFactorGraph& graph) const;

    // ISAM2Result firstupdatestep(NonlinearFactorGraph& newFactors,std::vector<int>& removeFactorIndices,
    //                             NonlinearFactorGraph& removeFactors,std::set<int>& unusedKeys,
    //                         std::set<int>& unusedIndices,  bool force_relinearize,bool* relinearizestep);
//functions from Bayestree;

protected:

    std::set<int> getAffectedFactors(const std::list<int>& keys,const ISAM2Data& isam2data) const;
    GaussianFactorGraph* relinearizeAffectedFactors(const std::list<int>& affectedKeys,
    const std::set<int>& relinKeys,const ISAM2Data& isam2data,std::set<int>& cachedlinearized) const;
//    GaussianFactorGraph* getCachedBoundaryFactors(std::list<ISAM2CliquePointer*>* orphans);
    int getCachedBoundaryFactors(std::list<ISAM2CliquePointer*>* orphans,GaussianFactorGraph& cachedBoundary);

    std::set<int> recalculate(NonlinearFactorGraph& newFactors,const std::set<int>& markedKeys, const std::set<int>& relinKeys,
                              const std::vector<int>& observedKeys, const std::set<int>& unusedIndices,
                              const std::map<int,int> *constrainKeys,
                              ISAM2Result& result,ISAM2Data& isam2data);
    void updateDelta(ISAM2Data& isam2data, bool forceFullSolve = false);

protected: //functions from BayesTree

    /** private helper method for saving the Tree to a text file in GraphViz format */
    void saveGraph(std::ostream &s,const ISAM2CliquePointer* clique,
                   int parentnum = 0) const;

    /** Gather data on a single clique */
    //  void getCliqueData(BayesTreeCliqueData stats, const ISAM2Clique& clique) const;

    /** remove a clique: warning, can result in a forest */
    int removeClique(ISAM2CliquePointer* clique);

    /** Fill the nodes index for a subtree */
    //void fillNodesIndex(const ISAM2Clique& subtree);

    // Friend JunctionTree because it directly fills roots and nodes index.
    friend class EliminatableClusterTree;

}; // ISAM2

#endif // ISAM2POINTER_H_INCLUDED
