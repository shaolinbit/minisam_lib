#ifndef ORDERING_H
#define ORDERING_H



/**
Metis order was discarded.
In order to improve the speed, Order class was deleted. Only the way to get Order is kept.
*/
#include "../inference/VariableIndex.h"
#include "../inference/FactorGraph.h"

#include <algorithm>
#include <vector>
namespace minisam
{

class VariableIndex;
class NonlinearFactorGraph;
 enum Ordering_OrderingType
    {
        Ordering_COLAMD,
        //METIS,
        Ordering_NATURAL, Ordering_CUSTOM
    };


 /// Invert (not reverse) the ordering - returns a map from key to order position
    std::map<int,int> Ordering_invert(const std::vector<int>& Ordering) ;

 /// Compute a fill-reducing ordering using COLAMD from a factor graph (see details for note on
    /// performance).
    template<class TPFactor>
    std::vector<int> Ordering_Colamd(const FactorGraph<TPFactor>& graph);
    /// Compute a fill-reducing ordering using COLAMD from a VariableIndex.
    std::vector<int> Ordering_Colamd(const VariableIndex& variableIndex);
    /// Compute a fill-reducing ordering using constrained COLAMD from a factor graph (see details
    /// for note on performance).  This internally builds a VariableIndex so if you already have a
    /// VariableIndex, it is faster to use COLAMD(const VariableIndex&).  This function constrains
    /// the variables in \c constrainLast to the end of the ordering, and orders all other variables
    /// before in a fill-reducing ordering.  If \c forceOrder is true, the variables in \c
    /// constrainLast will be ordered in the same order specified in the vector<Key> \c
    /// constrainLast.   If \c forceOrder is false, the variables in \c constrainLast will be
    /// ordered after all the others, but will be rearranged by CCOLAMD to reduce fill-in as well.

    template<class TPFactor>
    std::vector<int> Ordering_ColamdConstrainedLast(const FactorGraph<TPFactor>& graph,
                                   const std::vector<int>& constrainLast, bool forceOrder = false);
   /// Compute a fill-reducing ordering using constrained COLAMD from a VariableIndex.  This
    /// function constrains the variables in \c constrainLast to the end of the ordering, and orders
    /// all other variables before in a fill-reducing ordering.  If \c forceOrder is true, the
    /// variables in \c constrainLast will be ordered in the same order specified in the vector<Key>
    /// \c constrainLast.   If \c forceOrder is false, the variables in \c constrainLast will be
    /// ordered after all the others, but will be rearranged by CCOLAMD to reduce fill-in as well.
    std::vector<int> Ordering_ColamdConstrainedLast(
        const VariableIndex& variableIndex, const std::vector<int>& constrainLast,
        bool forceOrder = false);

   /// Compute a fill-reducing ordering using constrained COLAMD from a factor graph (see details
    /// for note on performance).  This internally builds a VariableIndex so if you already have a
    /// VariableIndex, it is faster to use COLAMD(const VariableIndex&).  In this function, a group
    /// for each variable should be specified in \c groups, and each group of variables will appear
    /// in the ordering in group index order.  \c groups should be a map from Key to group index.
    /// The group indices used should be consecutive starting at 0, but may appear in \c groups in
    /// arbitrary order.  Any variables not present in \c groups will be assigned to group 0.  This
    /// function simply fills the \c cmember argument to CCOLAMD with the supplied indices, see the
    /// CCOLAMD documentation for more information.

    template<class TPFactor>
    std::vector<int> Ordering_ColamdConstrained(const FactorGraph<TPFactor>& graph,
                               const std::map<int, int>& groups);
      /// Compute a fill-reducing ordering using constrained COLAMD from a VariableIndex.  In this
    /// function, a group for each variable should be specified in \c groups, and each group of
    /// variables will appear in the ordering in group index order.  \c groups should be a map from
    /// Key to group index. The group indices used should be consecutive starting at 0, but may
    /// appear in \c groups in arbitrary order.  Any variables not present in \c groups will be
    /// assigned to group 0.  This function simply fills the \c cmember argument to CCOLAMD with the
    /// supplied indices, see the CCOLAMD documentation for more information.
    std::vector<int> Ordering_ColamdConstrained(
        const VariableIndex& variableIndex, const std::map<int, int>& groups);
    std::vector<int> Ordering_ColamdConstrained(const VariableIndex& variableIndex,
                                     std::vector<int>& cmember);

    template<class TPFactor>
    std::vector<int> Ordering_Natural(const FactorGraph<TPFactor> &fg);

    template<class TPFactor>
     std::vector<int> Ordering_Create(Ordering_OrderingType orderingType,
                    const FactorGraph<TPFactor>& graph);
     std::vector<int> Ordering_Create(Ordering_OrderingType orderingType,
                    const NonlinearFactorGraph& graph);
        /// Compute a fill-reducing ordering using constrained COLAMD from a factor graph (see details
    /// for note on performance).  This internally builds a VariableIndex so if you already have a
    /// VariableIndex, it is faster to use COLAMD(const VariableIndex&).  This function constrains
    /// the variables in \c constrainLast to the end of the ordering, and orders all other variables
    /// before in a fill-reducing ordering.  If \c forceOrder is true, the variables in \c
    /// constrainFirst will be ordered in the same order specified in the vector<Key> \c
    /// constrainFirst.   If \c forceOrder is false, the variables in \c constrainFirst will be
    /// ordered before all the others, but will be rearranged by CCOLAMD to reduce fill-in as well.
    template<class TPFactor>
    std::vector<int> Ordering_ColamdConstrainedFirst(const FactorGraph<TPFactor>& graph,
                                    const std::vector<int>& constrainFirst, bool forceOrder = false);


    /// Compute a fill-reducing ordering using constrained COLAMD from a VariableIndex.  This
    /// function constrains the variables in \c constrainFirst to the front of the ordering, and
    /// orders all other variables after in a fill-reducing ordering.  If \c forceOrder is true, the
    /// variables in \c constrainFirst will be ordered in the same order specified in the
    /// vector<Key> \c constrainFirst.   If \c forceOrder is false, the variables in \c
    /// constrainFirst will be ordered before all the others, but will be rearranged by CCOLAMD to
    /// reduce fill-in as well.
    std::vector<int> Ordering_ColamdConstrainedFirst(
        const VariableIndex& variableIndex,
        const std::vector<int>& constrainFirst, bool forceOrder = false);

};
#endif // ORDERING_H

