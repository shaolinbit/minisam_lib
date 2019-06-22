#ifndef ISAM2JUNCTIONTREE_H_INCLUDED
#define ISAM2JUNCTIONTREE_H_INCLUDED


//EliminatableClusterTree support isam2.


#include "../inference/ClusterTree.h"
#include "../inference/EliminationTree.h"
#include "../symbolic/SymbolicConditional.h"

namespace minisam
{
class EliminationTree;
/**
 ISAM2JunctionTree is a kind of junctiontree using isam2 as its bayestree structure.
 */

class ISAM2JunctionTree : public EliminatableClusterTree
{
public:

    /// @name Standard Constructors
    /// @{
    /** Build the junction tree from an elimination tree. */
    ISAM2JunctionTree(const EliminationTree& eliminationTree,const GaussianFactorGraph& gf);

    /// @}

private:

    // Private default constructor (used in static construction methods)
    ISAM2JunctionTree() {}

};
};

#endif // ISAM2JUNCTIONTREE_H_INCLUDED
