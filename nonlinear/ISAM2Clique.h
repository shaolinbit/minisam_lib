#ifndef ISAM2CLIQUE_H_INCLUDED
#define ISAM2CLIQUE_H_INCLUDED

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    ISAM2Clique.h
 * @brief   Specialized iSAM2 Clique
 * @author  Michael Kaess, Richard Roberts
 */


#include "../linear/GaussianFactorGraph.h"
#include "../linear/GaussianConditional.h"
#include "../inference/EliminationTree.h"
#include "../inference/Symbol.h"
namespace minisam
{
/**
 * Specialized Clique structure for ISAM2, incorporating caching and gradient contribution
 * TODO: more documentation
 */
class ISAM2Clique
{
public:
    GaussianConditional* conditional_;
    ISAM2Clique* parent_;
    std::vector<ISAM2Clique*>* children_;
    GaussianFactorGraph* cachedSeparatorMarginal_;
    std::map<int,int>* assignedkeys;
    int problemSize_;

    RealGaussianFactor* cachedFactor_;
    std::map<int, std::map<int,Eigen::VectorXd>::iterator>* solnPointers_;

    bool setErased;
    bool isroot;

    /// Default constructor
    ISAM2Clique();// : Base() {}
    ~ISAM2Clique();


    /// Copy constructor, does *not* copy solution pointers as these are invalid in different trees.
    ISAM2Clique(const ISAM2Clique& other) :
        conditional_(other.conditional_),
        assignedkeys(other.assignedkeys),
        problemSize_(other.problemSize_),
        children_(other.children_),
        cachedSeparatorMarginal_(other.cachedSeparatorMarginal_), cachedFactor_(other.cachedFactor_),//  gradientContribution_(other.gradientContribution_),
        solnPointers_(other.solnPointers_),
        setErased(other.setErased) {}

    /// Assignment operator, does *not* copy solution pointers as these are invalid in different trees.
    ISAM2Clique operator=(const ISAM2Clique& other);

    //*//from BTC
    /** Access the conditional */
    GaussianConditional* conditional() const;

    /** is this the root of a Bayes tree ? */
    inline bool isRoot() const
    {
        return isroot;
    }
    bool operator==(const ISAM2Clique& other)const;


    /** Problem size (used for parallel traversal) */
    int problemSize() const;


    /** return a shared_ptr to the parent clique */
    ISAM2Clique* parent() const;

    GaussianFactorGraph* cachedSeparatorMarginal();
    //*//from BTC

    /// Overridden to also store the remaining factor and gradient contribution
    void setEliminationResult(const std::pair<GaussianConditional*, RealGaussianFactor*>& eliminationResult);

    /** Access the cached factor */
    RealGaussianFactor* cachedFactor();

    /// @}
    /// @name Advanced Interface
    /// @{

    /**
     * This deletes the cached shortcuts of all cliques (subtree) below this clique.
     * This is performed when the bayes tree is modified.
     */
    void deleteCachedShortcuts();

    ///@}

protected:

    /// Calculate set \f$ S \setminus B \f$ for shortcut calculations
    std::vector<int> separator_setminus_B(const ISAM2Clique* B) const;


};

class ISAM2OrphanWrapper: public JacobianFactor
{
public:
    //
    ISAM2Clique* clique_;
public:
    ISAM2OrphanWrapper(ISAM2Clique* clique,const std::vector<int>& inkeys);
};
};

#endif // ISAM2CLIQUE_H_INCLUDED
