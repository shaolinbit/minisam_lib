#ifndef ISAM2CLIQUEPOINTER_H_INCLUDED
#define ISAM2CLIQUEPOINTER_H_INCLUDED


#include "../linear/GaussianFactorGraph.h"
#include "../linear/GaussianConditional.h"
#include "../linear/GaussianBayesNetPointer.h"
#include "../inference/EliminationTree.h"
#include "../inference/Symbol.h"

/**
 * Specialized Clique structure for ISAM2, incorporating caching and gradient contribution
 * TODO: more documentation
 */
class ISAM2CliquePointer
{
public:
    GaussianConditional* conditional_;
    //int parentindex_;
    ISAM2CliquePointer* parent_;
    //int currentindex_;
    //std::vector<int> childrenindex;
    std::vector<ISAM2CliquePointer*>* children_;
    GaussianFactorGraph* cachedSeparatorMarginal_;
    std::map<int,int>* assignedkeys;
    int problemSize_;

    RealGaussianFactor* cachedFactor_;
   // Eigen::VectorXd* gradientContribution_;
    std::map<int, std::map<int,Eigen::VectorXd>::iterator>* solnPointers_;

    bool setErased;
    bool isroot;

    /// Default constructor
    ISAM2CliquePointer();// : Base() {}
    ~ISAM2CliquePointer();

    /*ISAM2CliquePointer(int parentindex,int currentindex)
    {
        parentindex_=parentindex;
        currentindex_=currentindex;
         conditional_=NULL;
    cachedSeparatorMarginal_=NULL;
    cachedFactor_=NULL;
        setErased=false;
    }*/

    /// Copy constructor, does *not* copy solution pointers as these are invalid in different trees.
    ISAM2CliquePointer(const ISAM2CliquePointer& other) :
        conditional_(other.conditional_),
        assignedkeys(other.assignedkeys),
        problemSize_(other.problemSize_),
        children_(other.children_),
        cachedSeparatorMarginal_(other.cachedSeparatorMarginal_), cachedFactor_(other.cachedFactor_),//  gradientContribution_(other.gradientContribution_),
      solnPointers_(other.solnPointers_),
        setErased(other.setErased){}

    /// Assignment operator, does *not* copy solution pointers as these are invalid in different trees.
    ISAM2CliquePointer operator=(const ISAM2CliquePointer& other)
    {
        conditional_=other.conditional_;
        problemSize_=other.problemSize_;
        parent_=other.parent_;
        children_=other.children_;
        assignedkeys=other.assignedkeys;
        cachedSeparatorMarginal_=other.cachedSeparatorMarginal_;
        cachedFactor_ = other.cachedFactor_;
        //gradientContribution_ = other.gradientContribution_;
        solnPointers_=other.solnPointers_;
        setErased=other.setErased;
        return *this;
    }

    //*//from BTC
    /** Access the conditional */
    GaussianConditional* conditional() const
    {
        return conditional_;
    }

    /** is this the root of a Bayes tree ? */
    inline bool isRoot() const {
        return isroot;
         }

    bool operator==(const ISAM2CliquePointer& other)const
    {
        return (conditional_==other.conditional_);
    }


    /** Problem size (used for parallel traversal) */
    int problemSize() const
    {
        return problemSize_;
    }


    /** return a shared_ptr to the parent clique */
    ISAM2CliquePointer* parent() const
    {
        return parent_;
    }

    GaussianFactorGraph* cachedSeparatorMarginal()
    {
        return cachedSeparatorMarginal_;
    }
    //*//from BTC

    /// Overridden to also store the remaining factor and gradient contribution
    void setEliminationResult(const std::pair<GaussianConditional*, RealGaussianFactor*>& eliminationResult);

    /** Access the cached factor */
    RealGaussianFactor* cachedFactor()
    {
        return cachedFactor_;
    }

    /** Access the gradient contribution
    const Eigen::VectorXd gradientContribution() const
    {
        return *gradientContribution_;
    }*/

    /** The size of subtree rooted at this clique, i.e., nr of Cliques */
//   int treeSize() const;

    /** Collect number of cliques with cached separator marginals */
//    int numCachedSeparatorMarginals() const;


    /// @}
    /// @name Advanced Interface
    /// @{

    /** return the conditional P(S|Root) on the separator given the root */
    GaussianBayesNetPointer* shortcut(const ISAM2CliquePointer* root,
                                      int EliminationFunctionType) const;

    /** return the marginal P(S) on the separator */
    GaussianFactorGraph* separatorMarginal(int EliminationFunctionType);

    /** return the marginal P(C) of the clique, using marginal caching */
    GaussianFactorGraph* marginal2(const int EliminationFunctionType);

    /**
     * This deletes the cached shortcuts of all cliques (subtree) below this clique.
     * This is performed when the bayes tree is modified.
     */
    void deleteCachedShortcuts();



    //friend class BayesTree;

protected:

    /// Calculate set \f$ S \setminus B \f$ for shortcut calculations
    std::vector<int> separator_setminus_B(const ISAM2CliquePointer* B) const;

    /** Determine variable indices to keep in recursive separator shortcut calculation The factor
     *  graph p_Cp_B has keys from the parent clique Cp and from B. But we only keep the variables
     *  not in S union B. */
    std::vector<int> shortcut_indices(const ISAM2CliquePointer* B, const GaussianFactorGraph& p_Cp_B) const;

public:
    /** Non-recursive delete cached shortcuts and marginals - internal only. */
    void deleteCachedShortcutsNonRecursive()
    {
        cachedSeparatorMarginal_->factors_.clear();
    }

};

class ISAM2OrphanWrapperPointer : public JacobianFactor
{
public:
    //
     ISAM2CliquePointer* clique_;
public:
    ISAM2OrphanWrapperPointer(ISAM2CliquePointer* clique,const std::vector<int>& inkeys):JacobianFactor(JacobianFactor())
    {
        // Store parent keys in our base type factor so that eliminating those parent keys will pull
        // this subtree into the elimination.
        this->keys_=inkeys;
       // this->cliqueindex_=cliqueindex;
       clique_=clique;
       this->iswrapper_=true;
    }
};


#endif // ISAM2CLIQUEPOINTER_H_INCLUDED
