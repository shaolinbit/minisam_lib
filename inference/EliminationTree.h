#ifndef ELIMINATIONTREE_H_INCLUDED
#define ELIMINATIONTREE_H_INCLUDED

/**
* @file    EliminationTree.h
*/
#include <utility>
#include "../linear/GaussianFactorGraph.h"
#include "../linear/RealGaussianFactor.h"
#include "../linear/GaussianBayesNet.h"


namespace minisam
{
class VariableIndex;
class GaussianFactorGraph;

/**
* An elimination tree is a data structure used intermediately during
* elimination.  In future versions it will be used to save work between
* multiple eliminations.
*
* When a variable is eliminated, a new factor is created by combining that
* variable's neighboring factors.  The new combined factor involves the combined
* factors' involved variables.  When the lowest-ordered one of those variables
* is eliminated, it consumes that combined factor.  In the elimination tree,
* that lowest-ordered variable is the parent of the variable that was eliminated to
* produce the combined factor.  This yields a tree in general, and not a chain
* because of the implicit sparse structure of the resulting Bayes net.
*
* This structure is examined even more closely in a JunctionTree, which
* additionally identifies cliques in the chordal Bayes net.
* \nosubgrouping
*/
struct ETNode
{
    int key; ///< key associated with root
    std::vector<int>* factorsindex; ///< factors associated with root
    std::vector<int>* childrenindex;
    ETNode()
    {
        factorsindex=new std::vector<int>();
        childrenindex=new std::vector<int>();
    }
    ~ETNode()
    {
        if(factorsindex!=NULL)
        {
            factorsindex->resize(0);
            delete factorsindex;
            factorsindex=NULL;
        }
        if(childrenindex!=NULL)
        {
            childrenindex->resize(0);
            delete childrenindex;
            childrenindex=NULL;
        }
    }
    RealGaussianFactor* eliminate(GaussianBayesNet* output,
                                  const std::vector<RealGaussianFactor*>& childrenFactors,
                                  const GaussianFactorGraph& gf,
                                  const Factorization  Eliminatefunction=CHOLESKY);
};

class EliminationTree
{
public:
    std::vector<int>* roots_;
    std::vector<ETNode*>* nodes;
    std::vector<int>* remainingFactorsindex_;

    /// @name Standard Constructors
    /// @{

    /**
    * Build the elimination tree of a factor graph using pre-computed column structure.
    * @param factorGraph The factor graph for which to build the elimination tree
    * @param structure The set of factors involving each variable.  If this is not
    * precomputed, you can call the Create(const GaussianFactorGraph&)
    * named constructor instead.
    * @return The elimination tree
    */
public:
    EliminationTree(const GaussianFactorGraph& factorGraph,
                    const VariableIndex& structure, const std::vector<int>& order);
    /** Build the elimination tree of a factor graph.  Note that this has to compute the column
    * structure as a VariableIndex, so if you already have this precomputed, use the other
    * constructor instead.
    * @param factorGraph The factor graph for which to build the elimination tree
    */
    EliminationTree(const GaussianFactorGraph& factorGraph, const std::vector<int>& order);
    /** Copy constructor - makes a deep copy of the tree structure, but only pointers to factors are
     *  copied, factors are not cloned. */
    EliminationTree(const EliminationTree& other);

    /** Assignment operator - makes a deep copy of the tree structure, but only pointers to factors
     *  are copied, factors are not cloned. */
    EliminationTree& operator=(const EliminationTree& other);


    ~EliminationTree();

    /// @}

public:

    /** Eliminate the factors to a Bayes net and remaining factor graph
    * @param function The function to use to eliminate, see the namespace functions
    * in GaussianFactorGraph.h
    * @return The Bayes net and factor graph resulting from elimination
    */
    std::pair<GaussianBayesNet*, GaussianFactorGraph*>
    eliminate(const GaussianFactorGraph& gf,const Factorization Eliminatefunction=CHOLESKY) const;


    /// @name Testable
    /// @{

    /** Print the tree to cout */
    void print() const;

    ///@}

public:
    /// @name Advanced Interface
    /// @{
    /** Return the remaining factors that are not pulled into elimination */
    const std::vector<int> remainingFactorsindex() const;

    /** Swap the data of this tree with another one, this operation is very fast. */
    void swap(EliminationTree& other);

    ///@}

protected:
    /// Protected default constructor
    EliminationTree() {}

};
struct ETEliminationDataChildrenFactors
{
    int parentindex_;
    int currentindex_;
    std::vector<RealGaussianFactor*> childFactors;
    ETEliminationDataChildrenFactors(
        int pIndex,int cindex,int nChildren):
        parentindex_(pIndex),currentindex_(cindex)
    {
        childFactors.reserve(nChildren);
    }
    ~ETEliminationDataChildrenFactors()
    {
        for(RealGaussianFactor* bf:childFactors)
        {
            if(bf!=NULL)
            {
                if(bf->model_!=NULL)
                {
                  delete bf->model_;
                  bf->model_=NULL;
                }
                delete bf;
                bf=NULL;
            }

        }

    }
};

struct ETTraversalNode
{
    bool expanded;
    ETNode* treeNode;
    int parentindex;

    ETTraversalNode(ETNode* _treeNode, int _parentindex) :
        expanded(false), treeNode(_treeNode), parentindex(_parentindex)
    {
    }
};

std::vector<RealGaussianFactor*>
inferenceEliminateTree(GaussianBayesNet* result, const EliminationTree& tree,
                       const GaussianFactorGraph& gf,
                       const Factorization Eliminatefunction=CHOLESKY);
void ETDepthFirstForestBN(GaussianBayesNet* BNresult,
                          const EliminationTree& forest,
                          ETEliminationDataChildrenFactors* rootData,
                          const GaussianFactorGraph& gf,
                          const Factorization elimationtype=CHOLESKY);
};
#endif // ELIMINATIONTREE_H_INCLUDED
