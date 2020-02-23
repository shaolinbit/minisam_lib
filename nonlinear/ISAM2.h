#ifndef ISAM2POINTER_H_INCLUDED
#define ISAM2POINTER_H_INCLUDED

#include "../nonlinear/ISAM2-impl.h"
#include "../nonlinear/NonlinearFactorGraph.h"
#include "../nonlinear/ISAM2Params.h"
#include "../inference/VariableIndex.h"
#include "../nonlinear/ISAM2Clique.h"
#include "../nonlinear/ISAM2Data.h"


/**
1. Split the datastructure from isam2 class.
2. Only keep the minimum functions for isam2update.
3. Only support pose3 in this version.
*/

/**
 * @file    ISAM2.h
 * @brief   Incremental update functionality (ISAM2) for BayesTree, with fluid
 * relinearization.
 */
namespace minisam
{
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


class  ISAM2
{
public:

    std::map<int,ISAM2Clique*>* nodesbtc;
    std::vector<ISAM2Clique*>* roots_;
    /** The current parameters */
    ISAM2Params params_;
protected:

    int update_count_; ///< Counter incremented every update(), used to determine periodic relinearization
    /** The current Dogleg Delta (trust region radius) */
    mutable double doglegDelta_;
public:

    /** Create an empty ISAM2 instance */
    ISAM2(const ISAM2Params& params);

    /** Create an empty ISAM2 instance using the default set of parameters (see ISAM2Params) */
    ISAM2();

    /** default virtual destructor */
    ~ISAM2();


    ISAM2(const ISAM2& other);

    ISAM2& operator=(const ISAM2& other);

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
     * @return An ISAM2Result struct containing information about the update*/
    ISAM2Result update(NonlinearFactorGraph& newFactors,
                       const std::map<int,minimatrix*>& newTheta,
                       ISAM2Data& isam2data,
                       std::vector<int>* removeFactorIndices=NULL,
                       std::map<int,int> *constrainedKeys=NULL,
                       std::list<int> *noRelinKeys=NULL,
                       std::list<int> *extraReelimKeys=NULL,
                       bool force_relinearize = false);
    void calculateEstimate(ISAM2Data& isam2data);

    /** Compute an estimate for a single variable using its incomplete linear delta computed
     * during the last update.  This is faster than calling the no-argument version of
     * calculateEstimate, which operates on all variables.
     * @param key
     * @return
     */
    minimatrix* calculateEstimate(ISAM2Data& isam2data,int key);
    /** Access the current delta, computed during the last call to update */
    void getDelta(ISAM2Data& isam2data);

    std::map<int,minivector> gradientAtZero() const;

    /** Compute the linear error */
    double error(const std::map<int,minivector>& x) const;
    mutable int lastBacksubVariableCount;


    //functions from Bayestree;

    /** number of cliques */
    int size() const;

    /** Check if there are any cliques in the tree */
    inline bool empty() const
    {
        return (roots()->empty()&&nodesbtc->empty());
    }


    /** Access node by variable */
    ISAM2Clique* operator[](int j) const;
    ISAM2Clique* operator[](int j);

    /** return root cliques */
    std::vector<ISAM2Clique*>* roots() const;

    /** alternate syntax for matlab: find the clique that contains the variable with Key j */
    ISAM2Clique* clique(int j);

    /**
     * Read only with side effects
     */

    /** saves the Tree to a text file in GraphViz format */
    void saveGraph(std::ostream& stm) const;

    /// @}
    /// @name Advanced Interface
    /// @{

    /** Remove all nodes */
    void clear();
    void clearcliques();
    void clearallcliques();
    void clearall();

    /** Clear all shortcut caches - use before timing on marginal calculation to avoid residual cache data */
    void deleteCachedShortcuts();

    /**
     * Remove path from clique to root and return that path as factors
     * plus a list of orphaned subtree roots. Used in removeTop below.
     */
    void removePath(ISAM2Clique* clique, std::list<int>* affectedkeys,
                    std::list<ISAM2Clique*> *orphans);
    /**
     * Given a list of indices, turn "contaminated" part of the tree back into a factor graph.
     * Factors and orphans are added to the in/out arguments.
     */
    void removeTop(const std::vector<int>& keys, std::list<int>* affectedkeys, std::list<ISAM2Clique*> *orphans);

protected:

    std::set<int> getAffectedFactors(const std::list<int>& keys,const ISAM2Data& isam2data) const;
    GaussianFactorGraph* relinearizeAffectedFactors(const std::list<int>& affectedKeys,
            const std::set<int>& relinKeys, ISAM2Data& isam2data,std::set<int>& cachedlinearized) const;
    int getCachedBoundaryFactors(std::list<ISAM2Clique*>* orphans,GaussianFactorGraph& cachedBoundary);

    std::set<int> recalculate(NonlinearFactorGraph& newFactors,const std::set<int>& markedKeys, const std::set<int>& relinKeys,
                              const std::vector<int>& observedKeys, const std::set<int>& unusedIndices,
                              const std::map<int,int> *constrainKeys,
                              ISAM2Result& result,ISAM2Data& isam2data);
    void updateDelta(ISAM2Data& isam2data, bool forceFullSolve = false);

protected: //functions from BayesTree

    /** private helper method for saving the Tree to a text file in GraphViz format */
    void saveGraph(std::ostream &s,const ISAM2Clique* clique,
                   int parentnum = 0) const;

    /** remove a clique: warning, can result in a forest */
    int removeClique(ISAM2Clique* clique);

    // Friend JunctionTree because it directly fills roots and nodes index.
    friend class EliminatableClusterTree;

}; // ISAM2
};
#endif // ISAM2_H_INCLUDED
