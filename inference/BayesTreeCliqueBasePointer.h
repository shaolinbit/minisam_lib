#ifndef BAYESTREECLIQUEBASEPOINTER_H_INCLUDED
#define BAYESTREECLIQUEBASEPOINTER_H_INCLUDED


/**
 * @file    BayesTreeCliqueBase.h
 * @brief   Base class for cliques of a BayesTree
 * @author
 */

#include "../linear/GaussianFactorGraph.h"
#include "../linear/GaussianConditional.h"
#include "../linear/GaussianBayesNetPointer.h"
#include "../inference/EliminationTree.h"
#include "../inference/Symbol.h"

class BayesTreeCliqueBasePointer;
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


class BayesTreeCliqueBasePointer
{
public:
    BayesTreeCliqueBasePointer* parent_;
    GaussianConditional* conditional_;
    GaussianFactorGraph* cachedSeparatorMarginal_;
    std::vector<BayesTreeCliqueBasePointer*> children_;
    int problemSize_;
public:

    BayesTreeCliqueBasePointer();

    ~BayesTreeCliqueBasePointer();

    /** Construct from a conditional, leaving parent and child pointers uninitialized */
    BayesTreeCliqueBasePointer(GaussianConditional* conditional) :
        conditional_(conditional), problemSize_(1),
        parent_(NULL),
        cachedSeparatorMarginal_(new GaussianFactorGraph())
    {
        //setErased=false;
        children_.clear();
    }
    BayesTreeCliqueBasePointer(const BayesTreeCliqueBasePointer& btcb):conditional_(btcb.conditional_),
        problemSize_(btcb.problemSize_),parent_(btcb.parent_),
        children_(btcb.children_),//childrenindexpair(btcb.childrenindexpair),
        cachedSeparatorMarginal_(btcb.cachedSeparatorMarginal_)
    {}

    BayesTreeCliqueBasePointer operator=(const BayesTreeCliqueBasePointer& btcb)
    {
        conditional_=btcb.conditional_;
        problemSize_=btcb.problemSize_;
        parent_=btcb.parent_;
        children_=btcb.children_;
       // currentindex_=btcb.currentindex_;
      //  childrenindex=btcb.childrenindex;
        //childrenindexpair=btcb.childrenindexpair;
        // btcsym=btcb.btcsym;
        cachedSeparatorMarginal_=btcb.cachedSeparatorMarginal_;
        //assignedkeys=btcb.assignedkeys;
        //setErased=btcb.setErased;
        return *this;
    }


    //sharedConditional conditional_;
    bool operator==(const BayesTreeCliqueBasePointer& other)const
    {
        //return (this->btcsym.c_==other.btcsym.c_)&&(this->btcsym.j_==other.btcsym.j_);
        return (conditional_==other.conditional_);
    }


    /// Fill the elimination result produced during elimination.  Here this just stores the
    /// conditional and ignores the remaining factor, but this is overridden in ISAM2Clique
    /// to also cache the remaining factor.
    // void setEliminationResult(const typename FactorGraphType::EliminationResult& eliminationResult);
    void setEliminationResult(const std::pair<GaussianConditional*, RealGaussianFactor*>& eliminationResult);

    /** is this the root of a Bayes tree ? */

    inline bool isRoot() const {
         return parent_==NULL;
          }




    /** Problem size (used for parallel traversal) */
    int problemSize() const
    {
        return problemSize_;
    }


    /** return a ptr to the parent clique */
    BayesTreeCliqueBasePointer* parent() const
    {
        return parent_;
    }

    GaussianFactorGraph* cachedSeparatorMarginal()
    {
        return cachedSeparatorMarginal_;
    }


    /** The size of subtree rooted at this clique, i.e., nr of Cliques */
    // int treeSize() const;

    /** Collect number of cliques with cached separator marginals */
    //  int numCachedSeparatorMarginals() const;
    /// @}
    /// @name Advanced Interface
    /// @{

    /** return the conditional P(S|Root) on the separator given the root */
    // BayesNetType shortcut(const derived_ptr& root, Eliminate function = EliminationTraitsType::DefaultEliminate) const;
    GaussianBayesNetPointer* shortcut(const BayesTreeCliqueBasePointer* root,
                                      int EliminationFunctionType) const;

    /** return the marginal P(S) on the separator */
    GaussianFactorGraph* separatorMarginal(int EliminationFunctionType);

    /** return the marginal P(C) of the clique, using marginal caching */
    GaussianFactorGraph* marginal2(int EliminationFunctionType);

    /**
     * This deletes the cached shortcuts of all cliques (subtree) below this clique.
     * This is performed when the bayes tree is modified.
     */
    void deleteCachedShortcuts();

  int numCachedSeparatorMarginals() const;

    friend class BayesTreePointer;

protected:

    /// Calculate set \f$ S \setminus B \f$ for shortcut calculations
    std::vector<int> separator_setminus_B(const BayesTreeCliqueBasePointer* B) const;

    /** Determine variable indices to keep in recursive separator shortcut calculation The factor
     *  graph p_Cp_B has keys from the parent clique Cp and from B. But we only keep the variables
     *  not in S union B. */
    std::vector<int> shortcut_indices(const BayesTreeCliqueBasePointer* B, GaussianFactorGraph* p_Cp_B) const;

    /** Non-recursive delete cached shortcuts and marginals - internal only. */
    void deleteCachedShortcutsNonRecursive()
    {
        cachedSeparatorMarginal_->factors_.clear();
    }


};

//bool BTCBisEqual(BayesTreeCliqueBasePointer* b1,BayesTreeCliqueBasePointer* b2);
#endif // BAYESTREECLIQUEBASEPOINTER_H_INCLUDED
