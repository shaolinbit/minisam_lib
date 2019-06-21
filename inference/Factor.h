#ifndef FACTOR_H_INCLUDED
#define FACTOR_H_INCLUDED


/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    Factor.h
 * @brief   The base class for all factors
 * @author  Kai Ni
 * @author  Frank Dellaert
 * @author  Richard Roberts
 */


#include <vector>
#include <algorithm>


namespace minisam
{

/**
 * This is the base class for all factor types.  It is templated on a KEY type,
 * which will be the type used to label variables.  Key types currently in use
 * are Index with symbolic (IndexFactor, SymbolicFactorGraph) and
 * Gaussian factors (GaussianFactor, JacobianFactor, HessianFactor, GaussianFactorGraph),
 * and Key with nonlinear factors (NonlinearFactor, NonlinearFactorGraph).
 * though currently only IndexFactor and IndexConditional derive from this
 * class, using Index keys.  This class does not store any data other than its
 * keys.  Derived classes store data such as matrices and probability tables.
 *
 *
 * This class is \b not virtual for performance reasons - derived symbolic classes,
 * IndexFactor and IndexConditional, need to be created and destroyed quickly
 * during symbolic elimination.  GaussianFactor and NonlinearFactor are virtual.
 * \nosubgrouping
 */
class  Factor
{

public:
    /// Iterator over keys
    typedef std::vector<int>::iterator iterator;
    /// Const iterator over keys
    typedef std::vector<int>::const_iterator const_iterator;

    /// The keys involved in this factor
    std::vector<int> keys_;

    /// @name Standard Constructors
    /// @{

    /** Default constructor for I/O */
     Factor();
    virtual ~Factor();

    /** Construct factor from container of keys.  This constructor is used internally from derived factor
    *  constructors, either from a container of keys or from a boost::assign::list_of. */
     explicit Factor(const std::vector<int> keys) : keys_(keys.begin(), keys.end()) {}

     Factor(const int& key);

     Factor(const Factor& fc);
    /** Construct factor from iterator keys.  This constructor may be used internally from derived
    *  factor constructors, although our code currently does not use this. */
     Factor(std::vector<int>::iterator first, std::vector<int>::iterator last);
     Factor(std::vector<int>::const_iterator first, std::vector<int>::const_iterator last);
    /** Construct factor from container of keys.  This is called internally from derived factor static
    *  factor methods, as a workaround for not being able to call the protected constructors above. */
    static Factor FromKeys(const std::vector<int>& keys);

    /** Construct factor from iterator keys.  This is called internally from derived factor static
    *  factor methods, as a workaround for not being able to call the protected constructors above. */
    static Factor FromIterators(std::vector<int>::iterator first, std::vector<int>::iterator last);
//
    /// @}

public:
    /// @name Standard Interface
    /// @{

    /// First key
    int front() const;

    /// Last key
    int back() const;

    /// find
    std::vector<int>::const_iterator find(int key) const;

    /// Access the factor's involved variable keys
    const std::vector<int>& keys() const;
    /** Iterator at beginning of involved variable keys */
    std::vector<int>::const_iterator begin() const;

    /** Iterator at end of involved variable keys */
    std::vector<int>::const_iterator end() const;
    /**
    * @return the number of variables involved in this factor
    */
    int size() const;

    void push_back(int key);
    void clear();

    Factor& operator=(const Factor& rObj);

    /// @}


    /// @name Testable
    /// @{

protected:
    /// check equality
    bool equals(const Factor& other, double tol = 1e-9) const;

    /// @}

public:
    /// @name Advanced Interface
    /// @{

    /** @return keys involved in this factor */
    std::vector<int>& keys();
    /** Iterator at beginning of involved variable keys */
    std::vector<int>::iterator begin();

    /** Iterator at end of involved variable keys */
    std::vector<int>::iterator end();

    void erase(std::vector<int>::iterator needtoerase);
    void erase(int key);
    bool empty();

    bool operator==(const Factor& other)const;
    ///@}

};

};


#endif // FACTOR_H_INCLUDED
