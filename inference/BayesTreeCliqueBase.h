#ifndef BAYESTREECLIQUEBASE_H_INCLUDED
#define BAYESTREECLIQUEBASE_H_INCLUDED


/**
 * @file    BayesTreeCliqueBase.h
 * @brief   Base class for cliques of a BayesTree
 */

#include "../linear/GaussianFactorGraph.h"
#include "../linear/GaussianConditional.h"
#include "../linear/GaussianBayesNet.h"
#include "../inference/EliminationTree.h"
#include "../inference/Symbol.h"
namespace minisam
{

class BayesTreeCliqueBase;
/**
 * This is the base class for BayesTree cliques.  The default and standard derived type is
 * BayesTreeClique, but some algorithms, like iSAM2, use a different clique type in order to store
 * extra data along with the clique.
 *
 * This class is templated on the derived class (i.e. the curiously recursive template pattern).
 * The advantage of this over using virtual classes is that it avoids the need for casting to get
 * the derived type.  This is possible because all cliques in a BayesTree are the same type - if
 * they were not then we'd need a virtual class.
 *
 * @tparam DERIVED The derived clique type.
 * @tparam CONDITIONAL The conditional type.
 * \nosubgrouping */

class BayesTreeCliqueBase;
/**
 * This is the base class for BayesTree cliques.  The default and standard derived type is
 * BayesTreeClique, but some algorithms, like iSAM2, use a different clique type in order to store
 * extra data along with the clique.
 *
 * This class is templated on the derived class (i.e. the curiously recursive template pattern).
 * The advantage of this over using virtual classes is that it avoids the need for casting to get
 * the derived type.  This is possible because all cliques in a BayesTree are the same type - if
 * they were not then we'd need a virtual class.
 *
 * @tparam DERIVED The derived clique type.
 * @tparam CONDITIONAL The conditional type.
 * \nosubgrouping */


class BayesTreeCliqueBase
{
public:
    BayesTreeCliqueBase* parent_;
    GaussianConditional* conditional_;
    GaussianFactorGraph* cachedSeparatorMarginal_;
    std::vector<BayesTreeCliqueBase*>* children_;
    std::map<int,int>* assignedkeys;
    int problemSize_;
    bool setErased;
    bool isroot;

public:

    BayesTreeCliqueBase();

    ~BayesTreeCliqueBase();

    /** Construct from a conditional, leaving parent and child pointers uninitialized */
    BayesTreeCliqueBase(GaussianConditional* conditional);
    BayesTreeCliqueBase(const BayesTreeCliqueBase& btcb);

    BayesTreeCliqueBase operator=(const BayesTreeCliqueBase& btcb)
    {
        conditional_=btcb.conditional_;
        problemSize_=btcb.problemSize_;
        parent_=btcb.parent_;
        children_=btcb.children_;
        cachedSeparatorMarginal_=btcb.cachedSeparatorMarginal_;
        assignedkeys=btcb.assignedkeys;
        setErased=btcb.setErased;
        isroot=btcb.isroot;
        return *this;
    }


    bool operator==(const BayesTreeCliqueBase& other)const
    {
        return (conditional_==other.conditional_);
    }


    /// Fill the elimination result produced during elimination.  Here this just stores the
    /// conditional and ignores the remaining factor, but this is overridden in ISAM2Clique
    /// to also cache the remaining factor.
    void setEliminationResult(const std::pair<GaussianConditional*, RealGaussianFactor*>& eliminationResult);

    /** is this the root of a Bayes tree ? */

    inline bool isRoot() const
    {
        return parent_==NULL;
    }




    /** Problem size (used for parallel traversal) */
    int problemSize() const;
    /** return a ptr to the parent clique */
    BayesTreeCliqueBase* parent() const;

    GaussianFactorGraph* cachedSeparatorMarginal();

    /// @}
    /// @name Advanced Interface
    /// @{

    /** return the conditional P(S|Root) on the separator given the root */
    GaussianBayesNet* shortcut(const BayesTreeCliqueBase* root,
                               Factorization EliminationFunctionType=CHOLESKY) const;

    /** return the marginal P(S) on the separator */
    GaussianFactorGraph* separatorMarginal(Factorization EliminationFunctionType=CHOLESKY);

    /** return the marginal P(C) of the clique, using marginal caching */
    GaussianFactorGraph* marginal2(Factorization EliminationFunctionType=CHOLESKY);

    /**
     * This deletes the cached shortcuts of all cliques (subtree) below this clique.
     * This is performed when the bayes tree is modified.
     */
    void deleteCachedShortcuts();

    int numCachedSeparatorMarginals() const;

    friend class BayesTree;

protected:

    /// Calculate set \f$ S \setminus B \f$ for shortcut calculations
    std::vector<int> separator_setminus_B(const BayesTreeCliqueBase* B) const;

    /** Determine variable indices to keep in recursive separator shortcut calculation The factor
     *  graph p_Cp_B has keys from the parent clique Cp and from B. But we only keep the variables
     *  not in S union B. */
    std::vector<int> shortcut_indices(const BayesTreeCliqueBase* B, GaussianFactorGraph* p_Cp_B) const;

    /** Non-recursive delete cached shortcuts and marginals - internal only. */
    void deleteCachedShortcutsNonRecursive();


};

};
#endif // BAYESTREECLIQUEBASE_H_INCLUDED
