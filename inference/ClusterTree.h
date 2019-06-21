#ifndef CLUSTERTREEPOINTER_H_INCLUDED
#define CLUSTERTREEPOINTER_H_INCLUDED


///This file was modified a lot from clustertree.h in gtsam.

/**
 * @file EliminatableClusterTree.h
 * @date Oct 8, 2013
 * @author Kai Ni
 * @author Richard Roberts
 * @author Frank Dellaert
 * @brief Collects factorgraph fragments defined on variable clusters, arranged in a tree
 */

#include "../inference/Ordering.h"
#include "../linear/GaussianFactorGraph.h"
#include "../nonlinear/ISAM2.h"
#include "../symbolic/SymbolicConditional.h"

namespace minisam
{

/// A Cluster is just a collection of factors
class Cluster
{
public:
    std::vector<int> childrenclusterindex;
    int clusterindex;
    Ordering orderedFrontalKeys;  ///< Frontal keys of this node
    std::vector<int> factorsindex;  ///< Factors associated with this node
    int problemSize_;

    Cluster() : problemSize_(0),clusterindex(0) {}

    virtual ~Cluster() {}

    const int& operator[](int i) const
    {
        return childrenclusterindex[i];
    }

    /// Construct from factors associated with a single key
    template <class CONTAINER>
    Cluster(int key, const CONTAINER& factorsToAdd);

    /// Add factors associated with a single key
    template <class CONTAINER>
    void addFactors(int key, const CONTAINER& factorsToAdd);

    void addFactors(int key, const std::vector<int>& factorsToAdd);

    /// Add a child cluster
    void addChild(const Cluster* cluster);

    int nrChildren() const;

    int nrFactors() const;

    int nrFrontals() const;

    int problemSize() const;
    /// Return a vector with nrFrontal keys for each child
    std::vector<int> nrFrontalsOfChildren(const std::vector<Cluster*>& ctlist) const;

    /// Merge in given cluster
    void merge(const Cluster* cluster);

    /// Merge all children for which bit is set into this node
    void mergeChildren(const std::vector<bool>& merge,const std::vector<Cluster*>& ctlist);
};

/**
 * A cluster-tree is associated with a factor graph and is defined as in Koller-Friedman:
 * each node k represents a subset \f$ C_k \sub X \f$, and the tree is family preserving, in that
 * each factor \f$ f_i \f$ is associated with a single cluster and \f$ scope(f_i) \sub C_k \f$.
 * \nosubgrouping
 */
class ClusterTree
{
public:
    std::vector<int> roots_;
    std::vector<Cluster*> ctlist;

    /// @name Standard Constructors
    /// @{

    /** Copy constructor - makes a deep copy of the tree structure, but only pointers to factors are
     *  copied, factors are not cloned. */
    ClusterTree(const ClusterTree& other)
    {
        *this = other;
    }

    ~ClusterTree()
    {
        for(std::vector<Cluster*>::iterator bbi=ctlist.begin(); bbi!=ctlist.end(); bbi++)
        {
            delete *bbi;
            *bbi=NULL;
        }
    }

    /// @}

public:

    /// Default constructor
    ClusterTree() {}

    /// @name Advanced Interface
    /// @{

    void addRoot(const int& cluster)
    {
        roots_.push_back(cluster);
    }

    void addChildrenAsRoots(const Cluster* cluster)
    {
        for(int cindex:cluster->childrenclusterindex)
        {
            this->addRoot(cindex);
        }
    }

    int nrRoots() const
    {
        return roots_.size();
    }

    /** Return the set of roots (one for a tree, multiple for a forest) */
    const std::vector<int>& roots() const
    {
        return roots_;
    }

    /// @}
    /// @name Details

    /// Assignment operator - makes a deep copy of the tree structure, but only pointers to factors
    /// are copied, factors are not cloned.
    ClusterTree& operator=(const ClusterTree& other);
    /// print this ClusterTree
    void print() const;

    /// @}
};

/**
 * A cluster-tree that eliminates to a Bayes tree.
 */
class EliminatableClusterTree : public ClusterTree
{
protected:
    std::vector<int> remainingFactorsindex_;

    /// @name Standard Constructors
    /// @{

    /** Copy constructor - makes a deep copy of the tree structure, but only pointers to factors are
     *  copied, factors are not cloned. */
    EliminatableClusterTree(const EliminatableClusterTree& other) : ClusterTree(other)
    {
        *this = other;
    }

    /// @}

public:
    /// @name Standard Interface
    /// @{

    /** Eliminate the factors to a Bayes tree and remaining factor graph
     * @param function The function to use to eliminate, see the namespace functions
     * in GaussianFactorGraph.h
     * @return The Bayes tree and factor graph resulting from elimination
     */
    std::pair<ISAM2*, GaussianFactorGraph*> eliminateISAM2(const
            int Eliminatetype,std::list<ISAM2Clique*>* orphancliques,const GaussianFactorGraph& gf);

    /// @}
    /// @name Advanced Interface
    /// @{

    /** Return the remaining factors that are not pulled into elimination */
    const std::vector<int>& remainingFactorsindex() const
    {
        return remainingFactorsindex_;
    }

    /// @}

protected:
    /// @name Details

    /// Assignment operator - makes a deep copy of the tree structure, but only pointers to factors
    /// are copied, factors are not cloned.
    EliminatableClusterTree operator=(const EliminatableClusterTree& other);

    /// Default constructor to be used in derived classes
    EliminatableClusterTree() {}
    ~EliminatableClusterTree() {}

    /// @}
};
//}

struct EliminationDataISAM2childFactors
{

    std::vector<RealGaussianFactor*> childFactors;
    int myIndexInParent;
    ISAM2Clique* bayesTreeNode;
    EliminationDataISAM2childFactors* parentdata_;
    EliminationDataISAM2childFactors(EliminationDataISAM2childFactors* parentData):
        parentdata_(parentData)
    {
        bayesTreeNode=new ISAM2Clique();
        bayesTreeNode->setErased=false;
        childFactors.clear();
        if (parentData!=NULL)
        {
            myIndexInParent = parentData->childFactors.size();
            parentData->childFactors.push_back(new RealGaussianFactor());
        }
        else
        {
            myIndexInParent = 0;
        }
        // Set up BayesTree parent and child pointers
        if (parentData!=NULL)
        {
            if (parentData->parentdata_!=NULL) // If our parent is not the dummy node
                bayesTreeNode->parent_ = parentData->bayesTreeNode;
            parentData->bayesTreeNode->children_->push_back(bayesTreeNode);
        }
    }
};

void InitEliminationDataISAM2childFactors(EliminationDataISAM2childFactors *currentdata,
                                                 EliminationDataISAM2childFactors *parentdata);


RealGaussianFactor* I2EliminationPostOrderVisitor(Cluster* node,
        EliminationDataISAM2childFactors* myData,
        ISAM2* result,int Eliminatetype, std::list<ISAM2Clique*>* orphans,const GaussianFactorGraph& gf);

struct CTTraversalNodeISAM2
{
    bool expanded;
    Cluster* treeNode;
    EliminationDataISAM2childFactors* parentdata_;
    std::list<EliminationDataISAM2childFactors>::iterator dataPointer;

    CTTraversalNodeISAM2(Cluster* _treeNode, EliminationDataISAM2childFactors* parentdata, ISAM2Clique* parentclique) :
        expanded(false), treeNode(_treeNode), parentdata_(parentdata)
    {
    }
};

int I2ClusterTreeDepthFirstForest(EliminatableClusterTree* forest,
        EliminationDataISAM2childFactors* rootData,
        ISAM2* result,int EFunction, std::list<ISAM2Clique*>* orphans,const GaussianFactorGraph& gf);

class ConstructorTraversalDataChildFactors
{
public:
    int currentindex_;
    int parentindex_;
    Cluster* myJTNode;
    std::vector<SymbolicConditional*> childSymbolicConditionals;
    std::vector<Factor*> childSymbolicFactors;

public:
    ConstructorTraversalDataChildFactors(int parentindex,int currentindex) :
        currentindex_(currentindex),parentindex_(parentindex)
    {
        myJTNode=new Cluster();
        myJTNode->clusterindex=currentindex;
    }
    ~ConstructorTraversalDataChildFactors()
    {

        for(std::vector<SymbolicConditional*>::iterator bbi=childSymbolicConditionals.begin();
                bbi!=childSymbolicConditionals.end(); bbi++)
        {
            if(*bbi!=NULL)
            {
                delete *bbi;
                *bbi=NULL;
            }
        }

    }
};
};

#endif // CLUSTERTREEPOINTER_H_INCLUDED
