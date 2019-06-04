/**
 * @file    Conditional.h
 * @brief   Base class for conditional densities
 * @author
 */
#ifndef CONDITIONAL_H
#define CONDITIONAL_H

// \callgraph
//#pragma once

#include "../inference/Factor.h"


/**
 * TODO: Update comments. The following comments are out of date!!!
 *
 * Base class for conditional densities, templated on KEY type.  This class
 * provides storage for the keys involved in a conditional, and iterators and
 * access to the frontal and separator keys.
 *
 * Derived classes *must* redefine the Factor and shared_ptr typedefs to refer
 * to the associated factor type and shared_ptr type of the derived class.  See
 * IndexConditional and GaussianConditional for examples.
 * \nosubgrouping
 */
//template<class FACTOR, class DERIVEDCONDITIONAL>
template<class TPFactor>
class Conditional
{
public:
    /** The first nrFrontal variables are frontal and the rest are parents. */
    TPFactor nrFrontals_;
    TPFactor nrParents_;
    int sizenrFrontals;
    int sizenrParents;

public:
    /// Typedef to this class
    //typedef Conditional<FACTOR,DERIVEDCONDITIONAL> This;

    //public:
    /** View of the frontal keys (call frontals()) */
    // typedef boost::iterator_range<typename FACTOR::const_iterator> Frontals;

    /** View of the separator keys (call parents()) */
    //  typedef boost::iterator_range<typename FACTOR::const_iterator> Parents;

// protected:
    /// @name Standard Constructors
    /// @{

    /** Empty Constructor to make serialization possible */
    // Conditional() : nrFrontals_(0) {}

    /** Constructor */
    //Conditional(int nrFrontals) : nrFrontals_(nrFrontals) {}
    //Conditional(int nrFrontals) : nrFrontals_(nrFrontals) {}
    Conditional(const TPFactor& nrFrontals,const TPFactor& nrParents ):
        nrFrontals_(nrFrontals),nrParents_(nrParents)
    {
        //TPFactor nrFrontals_(nrFrontals.keys_);
        // TPFactor nrParents_(nrParents.keys_);
        sizenrFrontals=nrFrontals.size();
        sizenrParents=nrParents.size();
    }

    Conditional(TPFactor nrFrontals)
    {
        TPFactor nrFrontals_(nrFrontals.keys_);
        //TPFactor nrParents_(nrParents.keys_);
        sizenrFrontals=nrFrontals.size();
        nrParents_.clear();
        sizenrParents=0;
        // sizenrParents=nrParents.size();
    }
    /// @}
    /// @name Testable
    /// @{

    /// @}

public:
    /// @name Standard Interface
    /// @{
    Conditional();//{}

    ~Conditional();

    //  Conditional<Factor>::~Conditional();

    int getsizenrFrontals() const;// { return sizenrFrontals; }

    // int nrParents() const { return (*asFactor()).size() - nrFrontals_; }
    int getsizenrParents() const;// { return sizenrParents; }

    int firstFrontalKey() const;
    /*
    {
        if(sizenrFrontals>0)
        {
            return nrFrontals_.front();
        }
        else
        {
             throw std::invalid_argument("Requested Conditional::firstFrontalKey from a conditional with zero frontal keys");
        }
    }*/

    int firstParentsKey() const;
    /*
    {
       if(sizenrParents>0)
       {
           return nrParents_.front();
       }
       else
       {
            throw std::invalid_argument("Requested Conditional::firstFrontalKey from a conditional with zero frontal keys");
       }
    }*/
    /*
    int firstFrontalKey() const {
      if(nrFrontals_ > 0)
        return (*asFactor()).front();
      else
        throw std::invalid_argument("Requested Conditional::firstFrontalKey from a conditional with zero frontal keys");
    }
    */
    /*
        boost::iterator_range<std::vector<int>::iterator>  frontals() const { return boost::make_iterator_range(beginFrontals(), endFrontals()); }

        boost::iterator_range<std::vector<int>::iterator> parents() const { return boost::make_iterator_range(beginParents(), endParents()); }
    */
    std::vector<int>::const_iterator cbeginFrontals() const;// { return nrFrontals_.begin(); }


    std::vector<int>::const_iterator cendFrontals() const;//{ return nrFrontals_.end(); }

    std::vector<int>::const_iterator cbeginParents() const;// { return nrParents_.begin(); }

    std::vector<int>::const_iterator cendParents() const;// { return nrParents_.end(); }

    /// @}
    /// @name Advanced Interface
    /// @{

    TPFactor& nrFrontals()  const;//{ return nrFrontals_; }
    TPFactor& nrParents() const;//{ return &nrParents_; }


    std::vector<int>::iterator beginFrontals() const;// { return nrFrontals_.begin(); }


    std::vector<int>::iterator endFrontals() const;// { return nrFrontals_.end(); }

    std::vector<int>::iterator beginParents() const;// { return nrParents_.begin(); }

    std::vector<int>::iterator endParents() const;// { return nrParents_.end(); }

    //private:
    // Cast to factor type (non-const) (casts down to derived conditional type, then up to factor type)
    // Factor* asFactor() { return static_cast<Factor&>(static_cast<Conditional&>(*this)); }

    // Cast to derived type (const) (casts down to derived conditional type, then up to factor type)
    // const Factor* asFactor() const { return static_cast<const Factor&>(static_cast<const Conditional&>(*this)); }
    //*/

};

#endif // CONDITIONAL_H
