#ifndef SCATTER_H
#define SCATTER_H

/**
 * @file    Scatter.h
 * @brief   Maps global variable indices to slot indices
 */

#include "../linear/GaussianFactorGraph.h"

namespace minisam
{

class GaussianFactorGraph;

/// One SlotEntry stores the slot index for a variable, as well its dim.
struct  SlotEntry
{
    int key;
    int dimension;
    SlotEntry(int _key, int _dimension) : key(_key), dimension(_dimension) {}
    std::string toString() const;
    friend bool operator<(const SlotEntry& p, const SlotEntry& q)
    {
        return p.key < q.key;
    }
    static bool Zero(const SlotEntry& p)
    {
        return p.dimension==0;
    }
};

/**
 * Scatter is an intermediate data structure used when building a HessianFactor
 * incrementally, to get the keys in the right order. In spirit, it is a map
 * from global variable indices to slot indices in the union of involved
 * variables. We also include the dimensionality of the variable.
 */
class Scatter : public std::vector<SlotEntry>
{
public:
    /// Default Constructor
    Scatter() {}

    // Construct from gaussian factor graph, with optional (partial or complete) ordering
    Scatter(const GaussianFactorGraph& gfg,
            const std::vector<int>& ordering);
    Scatter(const std::vector<const RealGaussianFactor*>&factors, const std::vector<int>& ordering);
    Scatter(const GaussianFactorGraph& gfg);
    /// Add a key/dim pair
    void add(int key, int dim);

private:

    /// Find the SlotEntry with the right key (linear time worst case)
    iterator find(int key);
};
};
#endif // SCATTER_H
