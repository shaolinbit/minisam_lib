#ifndef NONLINEARISAMPOINTER_H_INCLUDED
#define NONLINEARISAMPOINTER_H_INCLUDED

/**
 * @file NonlinearISAM.h
 * @date
 * @author
 */

//#pragma once

#include "../nonlinear/NonlinearFactorGraph.h"
#include "../inference/ISAMPointer.h"
#include "../gmfconfig.h"
#include "../geometry/Pose2.h"
#include "../geometry/Pose3.h"

/**
 * Wrapper class to manage ISAM in a nonlinear context
 */
class  NonlinearISAMPointer
{
protected:

    /** The internal iSAM object */
    ISAMPointer isam_;


    /** The original factors, used when relinearizing */
    NonlinearFactorGraph factors_;

    /** The reordering interval and counter */
    int reorderInterval_;
    int reorderCounter_;

    /** The elimination function */
    int eliminationFunction_;
public:
    /** The current linearization point */
    std::map<int,Eigen::VectorXd> linPoint_vector;

#ifdef GMF_Using_Pose3
    std::map<int,Pose3> linPoint_pose;
#else
    std::map<int,Pose2> linPoint_pose;
#endif
public:

    /// @name Standard Constructors
    /// @{

    /**
     * Periodically reorder and relinearize
     * @param reorderInterval is the number of updates between reorderings,
     *   0 never reorders (and is dangerous for memory consumption)
     *  1 (default) reorders every time, in worse case is batch every update
     *  typical values are 50 or 100
     */
    NonlinearISAMPointer(int reorderInterval = 1,
                         const int eliminationFunction = 0) :
        reorderInterval_(reorderInterval), reorderCounter_(0), eliminationFunction_(eliminationFunction) {}

    /// @}
    /// @name Standard Interface
    /// @{

    /** Return the current solution estimate */
#ifdef GMF_Using_Pose3
    std::map<int,Eigen::VectorXd> estimate(std::map<int,Pose3>* pose3lin);
#else
    std::map<int,Eigen::VectorXd> estimate(std::map<int,Pose2>* pose2lin);
#endif // GMF_Using_Pose3
    /** find the marginal covariance for a single variable */

    Eigen::MatrixXd marginalCovariance(int key);


    // access

    /** access the underlying bayes tree */
    const ISAMPointer& bayesTreePointer() const
    {
        return isam_;
    }

    /** Return the current linearization point */
    const std::map<int,Eigen::VectorXd>& getLinearizationPointVector() const
    {
        return linPoint_vector;
    }

#ifdef GMF_Using_Pose3
    const std::map<int,Pose3>& getLinearizationPointPose() const
    {
        return linPoint_pose;
    }
#else
    const std::map<int,Pose2>& getLinearizationPointPose() const
    {
        return linPoint_pose;
    }
#endif

    /** get underlying nonlinear graph */
    const NonlinearFactorGraph& getFactorsUnsafe() const
    {
        return factors_;
    }

    /** get counters */
    int reorderInterval() const
    {
        return reorderInterval_;    ///<TODO: comment
    }
    int reorderCounter() const
    {
        return reorderCounter_;    ///<TODO: comment
    }

    /** saves the Tree to a text file in GraphViz format */
    void saveGraph(std::ostream& stm) const;

    /// @}
    /// @name Advanced Interface
    /// @{

    /** Add new factors along with their initial linearization points */
    void update(const NonlinearFactorGraph& newFactors, const std::map<int,Eigen::VectorXd>& initialValue);


#ifdef GMF_Using_Pose3
    void update(const NonlinearFactorGraph& newFactors,
                const std::map<int,Eigen::VectorXd>& initialValue,
                const std::map<int,Pose3>& initialPose);
#else
    void update(const NonlinearFactorGraph& newFactors,
                const std::map<int,Eigen::VectorXd>& initialValue,
                const std::map<int,Pose2>& initialPose);
#endif // GMF_Using_Pose3
    /** Relinearization and reordering of variables */
    void reorder_relinearize();

    /// @}

};



#endif // NONLINEARISAMPOINTER_H_INCLUDED
