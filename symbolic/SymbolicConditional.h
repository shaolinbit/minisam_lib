#ifndef SYMBOLICCONDITONAL_H
#define SYMBOLICCONDITONAL_H


/**
 * @file    SymbolicConditional.h
 * @author
 * @date
 */


#include "../inference/Factor.h"
//#include "../inference/FactorConditional.h"


/**
 * SymbolicConditional is a conditional with keys but no probability
 * data, produced by symbolic elimination of SymbolicFactor.
 *
 * It is also a SymbolicFactor, and thus derives from it.  It
 * derives also from Conditional<This>, which is a generic interface
 * class for conditionals.
 * \nosubgrouping
 */
class  SymbolicConditional :
    public Factor//,public FactorConditional
{
public:
    Factor* nrFrontals_;
    Factor* nrParents_;
public:

    /// @name Standard Constructors
    /// @{

    /** Empty Constructor to make serialization possible */
    SymbolicConditional();// {}

    /** No parents */
    // SymbolicConditional(Factor& ff,Factor& fFactor) : Factor(ff), Conditional(fFactor) {}

    SymbolicConditional(const Factor& fFactor) : Factor(fFactor), nrFrontals_(new Factor(fFactor)),nrParents_(new Factor()) {}

    SymbolicConditional(const Factor& fFactor,int sizeFrontal);

    SymbolicConditional(const SymbolicConditional& other);


    /** Single parent
    SymbolicConditional(int j, int parent) : BaseFactor(j, parent), BaseConditional(1) {}


    // Two parents
    SymbolicConditional(int j, int parent1, int parent2) : BaseFactor(j, parent1, parent2), BaseConditional(1) {}

    // Three parents
    SymbolicConditional(int j, int parent1, int parent2, int parent3) : BaseFactor(j, parent1, parent2, parent3), BaseConditional(1) {}

    // Named constructor from an arbitrary number of keys and frontals
    template<typename ITERATOR>
    static SymbolicConditional FromIterators(ITERATOR firstKey, ITERATOR lastKey, size_t nrFrontals)
    {
      SymbolicConditional result;
      (BaseFactor&)result = BaseFactor::FromIterators(firstKey, lastKey);
      result.nrFrontals_ = nrFrontals;
      return result;
    }

    // Named constructor from an arbitrary number of keys and frontals
    template<typename ITERATOR>
    static SymbolicConditional::shared_ptr FromIteratorsShared(ITERATOR firstKey, ITERATOR lastKey, size_t nrFrontals)
    {
      SymbolicConditional::shared_ptr result = boost::make_shared<SymbolicConditional>();
      result->keys_.assign(firstKey, lastKey);
      result->nrFrontals_ = nrFrontals;
      return result;
    }

     //Named constructor from an arbitrary number of keys and frontals
    template<class CONTAINER>
    static SymbolicConditional FromKeys(const CONTAINER& keys, size_t nrFrontals) {
      return FromIterators(keys.begin(), keys.end(), nrFrontals);
    }

    //Named constructor from an arbitrary number of keys and frontals
    template<class CONTAINER>
    static SymbolicConditional::shared_ptr FromKeysShared(const CONTAINER& keys, size_t nrFrontals) {
      return FromIteratorsShared(keys.begin(), keys.end(), nrFrontals);
    }
    */
    ~SymbolicConditional();// {}
  int getsizenrFrontals() const // { return sizenrFrontals; }
    {
        return nrFrontals_->size();
    }

    // int nrParents() const { return (*asFactor()).size() - nrFrontals_; }
    int getsizenrParents() const // { return sizenrParents; }
    {
        return nrParents_->size();
    }

    int firstFrontalKey() const;

    int firstParentsKey() const;
    std::vector<int>::const_iterator cbeginFrontals()
    const; // { return nrFrontals_.begin(); }

    std::vector<int>::const_iterator cendFrontals()
    const; //{ return nrFrontals_.end(); }

    std::vector<int>::const_iterator cbeginParents()
    const; // { return nrParents_.begin(); }

    std::vector<int>::const_iterator cendParents()
    const; // { return nrParents_.end(); }

    /// @}
    /// @name Advanced Interface
    /// @{

    Factor* nrFrontals() const; //{ return nrFrontals_; }
    Factor* nrParents() const;  //{ return &nrParents_; }

    /// Copy this object as its actual derived type.
    // SymbolicConditional& clone() const { return *this; }

    /// @}
};
#endif // SYMBOLICCONDITONAL_H
