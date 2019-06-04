#ifndef ISAMPOINTER_H_INCLUDED
#define ISAMPOINTER_H_INCLUDED


#include "../inference/BayesTreePointer.h"

/**
 * @file    ISAM.h
 * @brief   Incremental update functionality (iSAM) for BayesTree.
 * @author
 */


/**
 * A Bayes tree with an update methods that implements the iSAM algorithm.
 * Given a set of new factors, it re-eliminates the invalidated part of the tree.
 * \nosubgrouping
 */
class ISAMPointer: public BayesTreePointer
{
public:
    /// @name Standard Constructors
    /// @{

    /** Create an empty Bayes Tree */
    ISAMPointer() {}
    ~ISAMPointer() {}

    /** Copy constructor */
    ISAMPointer(const BayesTreePointer& bayesTree) : BayesTreePointer(bayesTree) {}

    /// @}
    /// @name Advanced Interface Interface
    /// @{

    /**
     * update the Bayes tree with a set of new factors, typically derived from measurements
     * @param newFactors is a factor graph that contains the new factors
     * @param function an elimination routine
     */
    void update(const GaussianFactorGraph& newFactors, const int EliminatefunctionType);

    /** update_internal provides access to list of orphans for drawing purposes */
    void update_internal(const GaussianFactorGraph& newFactors, std::list<BayesTreeCliqueBasePointer*>* orphans,
                         const int Eliminatefunction);

    /// @}

};



#endif // ISAMPOINTER_H_INCLUDED
