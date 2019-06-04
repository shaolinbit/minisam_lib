#ifndef CLUSTERTREEPOINTER_H_INCLUDED
#define CLUSTERTREEPOINTER_H_INCLUDED

/**
 * @file EliminatableClusterTree.h
 * @date Oct 8, 2013
 * @author
 * @brief Collects factorgraph fragments defined on variable clusters, arranged in a tree
 */

#include "../inference/Ordering.h"
#include "../linear/GaussianFactorGraph.h"
#include "../inference/BayesTreePointer.h"
#include "../nonlinear/ISAM2.h"
#include "../symbolic/SymbolicConditional.h"

//#include "../inference/JunctionTree.h"

//struct ConstructorTraversalDataChildFactors;
/// A Cluster is just a collection of factors
// TODO(frank): re-factor JunctionTree so we can make members private
class Cluster
{
public:
    std::vector<int> childrenclusterindex;
    int clusterindex;
    Ordering orderedFrontalKeys;  ///< Frontal keys of this node
    // GaussianFactorGraph factors;  ///< Factors associated with this node
    std::vector<int> factorsindex;  ///< Factors associated with this node
    int problemSize_;

    Cluster() : problemSize_(0),clusterindex(0) {}

    virtual ~Cluster() {}

    const int& operator[](int i) const
    {
        // for(Cluster& clusterc:*ctlist)
        // {
        //  if(childrenclusterindex[i]==clusterc.clusterindex)
        //    {
        return childrenclusterindex[i];
        //     }
        //  }
        // return (children[i]);
    }

    /// Construct from factors associated with a single key
    template <class CONTAINER>
    Cluster(int key, const CONTAINER& factorsToAdd)
        : problemSize_(0)
    {
        addFactors(key, factorsToAdd);
    }

    /// Add factors associated with a single key
    template <class CONTAINER>
    void addFactors(int key, const CONTAINER& factorsToAdd)
    {
        orderedFrontalKeys.push_back(key);
        factorsindex.push_back(factorsToAdd);
        problemSize_ += factorsindex.size();
    }

    void addFactors(int key, const std::vector<int>& factorsToAdd)
    {
        orderedFrontalKeys.push_back(key);
        //for(const RealGaussianFactor& ib:factorsToAdd)
        // factors.push_back(ib);
        for(int ib:factorsToAdd)
            factorsindex.push_back(ib);
        problemSize_ += factorsindex.size();
    }

    /// Add a child cluster
    void addChild(const Cluster* cluster)
    {
        childrenclusterindex.push_back(cluster->clusterindex);
        problemSize_ = std::max(problemSize_, cluster->problemSize_);
    }

    int nrChildren() const
    {
        return childrenclusterindex.size();
    }

    int nrFactors() const
    {
        return factorsindex.size();
    }

    int nrFrontals() const
    {
        return orderedFrontalKeys.size();
    }

    int problemSize() const
    {
        return problemSize_;
    }



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
    //std::vector<int> clusterindexvector;

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
      for(std::vector<Cluster*>::iterator bbi=ctlist.begin();bbi!=ctlist.end();bbi++)
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
        //for (auto child : cluster->children)
        //std::vector<Cluster>::const_iterator child=cluster.children.begin();
        /*for(std::vector<Cluster>::const_iterator child=ctlist.begin();
                child!=ctlist.end(); child++)
        {
            for(int cindex:cluster.childrenclusterindex)
            {
                if(cindex==child->clusterindex)
                {
                    this->addRoot(*child);
                    break;
                }
            }
        }*/
        for(int cindex:cluster->childrenclusterindex)
        {
            // this->addRoot(ctlist.at(cindex));
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

    /*
    const Cluster& operator[](int i) const
    {
        return roots_[i];
    }*/

    /// @}

//public:
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
    //std::vector<RealGaussianFactor> remainingFactors_;
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
    std::pair<BayesTreePointer*, GaussianFactorGraph*> eliminate(const
            int Eliminatetype,std::list<BayesTreeCliqueBasePointer*>* orphans,const GaussianFactorGraph& gf);

    std::pair<ISAM2*, GaussianFactorGraph*> eliminateISAM2(const
            int Eliminatetype,std::list<ISAM2CliquePointer*>* orphancliques,const GaussianFactorGraph& gf);

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
/* ************************************************************************* */
// Elimination traversal data - stores a pointer to the parent data and collects
// the factors resulting from elimination of the children.  Also sets up BayesTree
// cliques with parent and child pointers.
struct EliminationDatachildFactorsPointer
{
    //int parentindex_;
    //int currentindex_;

    int myIndexInParent;
    std::vector<RealGaussianFactor*> childFactors;


    BayesTreeCliqueBasePointer* bayesTreeNode;
    EliminationDatachildFactorsPointer* parentdata_;

    EliminationDatachildFactorsPointer(EliminationDatachildFactorsPointer* parentData):
        parentdata_(parentData)
    {
        bayesTreeNode=new BayesTreeCliqueBasePointer();
       // bayesTreeNode->currentindex_=currentindex;
       // bayesTreeNode->parentindex_=parentindex;
       // bayesTreeNode->setErased=false;
        childFactors.clear();
       // if(parentindex_==-1)
       //     myIndexInParent=0;
         if (parentData!=NULL) {
      myIndexInParent = parentData->childFactors.size();
      parentData->childFactors.push_back(new RealGaussianFactor());
    } else {
      myIndexInParent = 0;
    }
    // Set up BayesTree parent and child pointers
    if (parentData!=NULL) {
      if (parentData->parentdata_!=NULL) // If our parent is not the dummy node
        bayesTreeNode->parent_ = parentData->bayesTreeNode;
      parentData->bayesTreeNode->children_.push_back(bayesTreeNode);
    }
    }
};

void InitEliminationDatachildFactorsPointer(EliminationDatachildFactorsPointer *currentdata,EliminationDatachildFactorsPointer *parentdata);

struct EliminationDataISAM2childFactorsPointer
{

    std::vector<RealGaussianFactor*> childFactors;
    //int parentindex_;
    //int currentindex_;
    int myIndexInParent;
   // int childindex_;

    //BayesTreeCliqueBase bayesTreeNode;
    ISAM2CliquePointer* bayesTreeNode;
    EliminationDataISAM2childFactorsPointer* parentdata_;
    EliminationDataISAM2childFactorsPointer(EliminationDataISAM2childFactorsPointer* parentData):
        parentdata_(parentData)
    {
       // childindex_=0;
        bayesTreeNode=new ISAM2CliquePointer();
     //   bayesTreeNode->currentindex_=currentindex;
      //  bayesTreeNode->parentindex_=parentindex;
        bayesTreeNode->setErased=false;
       // if(parentindex_==-1)
      //      myIndexInParent=0;
       childFactors.clear();
       // if(parentindex_==-1)
       //     myIndexInParent=0;
         if (parentData!=NULL) {
      myIndexInParent = parentData->childFactors.size();
      parentData->childFactors.push_back(new RealGaussianFactor());
    } else {
      myIndexInParent = 0;
    }
    // Set up BayesTree parent and child pointers
    if (parentData!=NULL) {
      if (parentData->parentdata_!=NULL) // If our parent is not the dummy node
        bayesTreeNode->parent_ = parentData->bayesTreeNode;
      parentData->bayesTreeNode->children_->push_back(bayesTreeNode);
    }
    }
};

void InitEliminationDataISAM2childFactorsPointer(EliminationDataISAM2childFactorsPointer *currentdata,EliminationDataISAM2childFactorsPointer *parentdata);

// Elimination pre-order visitor - creates the EliminationData structure for the visited node.
static EliminationDatachildFactorsPointer* EliminationPreOrderVisitorPointer(
    const Cluster& node,
    EliminationDatachildFactorsPointer* parentData)
{
    //assert(node);
    EliminationDatachildFactorsPointer* myData=new EliminationDatachildFactorsPointer(parentData);
    InitEliminationDatachildFactorsPointer(myData,parentData);
    myData->bayesTreeNode->problemSize_ = node.problemSize();
    return myData;
}

// Elimination pre-order visitor - creates the EliminationData structure for the visited node.
static EliminationDataISAM2childFactorsPointer* EliminationPreOrderVisitorPointer(
    const Cluster& node,
    EliminationDataISAM2childFactorsPointer* parentData)
{
    //assert(node);
    EliminationDataISAM2childFactorsPointer* myData=new EliminationDataISAM2childFactorsPointer(parentData);
    InitEliminationDataISAM2childFactorsPointer(myData,parentData);
    myData->bayesTreeNode->problemSize_ = node.problemSize();
    return myData;
}


// Elimination post-order visitor - combine the child factors with our own factors, add the
// resulting conditional to the BayesTree, and add the remaining factor to the parent.


// Elimination post-order visitor - combine the child factors with our own factors, add the
// resulting conditional to the BayesTree, and add the remaining factor to the parent.

RealGaussianFactor* EliminationPostOrderVisitorPointer(Cluster* node, EliminationDatachildFactorsPointer* myData,
        BayesTreePointer* result,int Eliminatetype, std::list<BayesTreeCliqueBasePointer*>* orphans,const GaussianFactorGraph& gf);
RealGaussianFactor* I2EliminationPostOrderVisitorPointer(Cluster* node,
        EliminationDataISAM2childFactorsPointer* myData,
        ISAM2* result,int Eliminatetype, std::list<ISAM2CliquePointer*>* orphans,const GaussianFactorGraph& gf);
struct CTTraversalNodePointer
{
    bool expanded;
    Cluster* treeNode;
   // int parentindex;
    EliminationDatachildFactorsPointer* parentdata_;
    std::list<EliminationDatachildFactorsPointer>::iterator dataPointer;
 //   BayesTreeCliqueBasePointer* parentclique_;


    CTTraversalNodePointer(Cluster* _treeNode, EliminationDatachildFactorsPointer* parentdata, BayesTreeCliqueBasePointer* parentclique) :
        expanded(false), treeNode(_treeNode), parentdata_(parentdata)//,parentclique_(parentclique)
    {
    }
};

struct CTTraversalNodeISAM2Pointer
{
    bool expanded;
    Cluster* treeNode;
    //int parentindex;
    EliminationDataISAM2childFactorsPointer* parentdata_;
    std::list<EliminationDataISAM2childFactorsPointer>::iterator dataPointer;
  //  ISAM2CliquePointer* parentclique_;

    CTTraversalNodeISAM2Pointer(Cluster* _treeNode, EliminationDataISAM2childFactorsPointer* parentdata, ISAM2CliquePointer* parentclique) :
        expanded(false), treeNode(_treeNode), parentdata_(parentdata)//,parentclique_(parentclique)
    {
    }
};


//void CTDepthFirstForestClone(const ClusterTree& forest, Cluster*  rootContainer);

int ClusterTreeDepthFirstForestPointer(EliminatableClusterTree* forest, EliminationDatachildFactorsPointer* rootData,
                                        BayesTreePointer* result,int EFunction, std::list<BayesTreeCliqueBasePointer*>* orphans,
                                        const GaussianFactorGraph& gf);

int I2ClusterTreeDepthFirstForestPointer(EliminatableClusterTree* forest,
        EliminationDataISAM2childFactorsPointer* rootData,
        ISAM2* result,int EFunction, std::list<ISAM2CliquePointer*>* orphans,const GaussianFactorGraph& gf);

class ConstructorTraversalDataChildFactorsPointer
{
public:
    int currentindex_;
    int parentindex_;
    Cluster* myJTNode;
    std::vector<SymbolicConditional*> childSymbolicConditionals;
    std::vector<Factor*> childSymbolicFactors;

public:
    ConstructorTraversalDataChildFactorsPointer(int parentindex,int currentindex) :
        currentindex_(currentindex),parentindex_(parentindex)
    {
        myJTNode=new Cluster();
        myJTNode->clusterindex=currentindex;
    }
    ~ConstructorTraversalDataChildFactorsPointer()
    {

        for(std::vector<SymbolicConditional*>::iterator bbi=childSymbolicConditionals.begin();
        bbi!=childSymbolicConditionals.end();bbi++)
        {
            if(*bbi!=NULL)
            {
            delete *bbi;
            *bbi=NULL;
            }
        }

    }
};


#endif // CLUSTERTREEPOINTER_H_INCLUDED
