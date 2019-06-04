#ifndef FACTOR_H_INCLUDED
#define FACTOR_H_INCLUDED

#include <vector>
#include <algorithm>

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
 * Note that derived classes *must* redefine the ConditionalType and shared_ptr
 * typedefs to refer to the associated conditional and shared_ptr types of the
 * derived class.  See IndexFactor, JacobianFactor, etc. for examples.
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
    Factor()
    {
        keys_.clear();
    }
    virtual ~Factor();

    /** Construct factor from container of keys.  This constructor is used internally from derived factor
    *  constructors, either from a container of keys or from a boost::assign::list_of. */
    //template<typename CONTAINER>
    explicit Factor(const std::vector<int> keys) : keys_(keys.begin(), keys.end()) {}

    Factor(const int& key)
    {
        keys_.clear();
        keys_.push_back(key);
    }

    Factor(const Factor& fc):keys_(fc.keys()) {}

    /** Construct factor from iterator keys.  This constructor may be used internally from derived
    *  factor constructors, although our code currently does not use this. */
    //  template<typename ITERATOR>
    Factor(std::vector<int>::iterator first, std::vector<int>::iterator last) : keys_(first, last) {}
    Factor(std::vector<int>::const_iterator first, std::vector<int>::const_iterator last) : keys_(first, last) {}
    /** Construct factor from container of keys.  This is called internally from derived factor static
    *  factor methods, as a workaround for not being able to call the protected constructors above. */
    //template<typename CONTAINER>
    static Factor FromKeys(const std::vector<int>& keys)
    {
        return Factor(keys);
    }

    /** Construct factor from iterator keys.  This is called internally from derived factor static
    *  factor methods, as a workaround for not being able to call the protected constructors above. */
    static Factor FromIterators(std::vector<int>::iterator first, std::vector<int>::iterator last)
    {
        return Factor(first, last);
    }
//
    /// @}

    //virtual Factor* clone() const=0;

public:
    /// @name Standard Interface
    /// @{

    /// First key
    int front() const
    {
        return keys_.front();
    }

    /// Last key
    int back() const
    {
        return keys_.back();
    }

    /// find
    std::vector<int>::const_iterator find(int key) const
    {
        return std::find(begin(), end(), key);
    }

    /// Access the factor's involved variable keys
    const std::vector<int>& keys() const
    {
        return keys_;
    }
    /** Iterator at beginning of involved variable keys */
    std::vector<int>::const_iterator begin() const
    {
        return keys_.begin();
    }

    /** Iterator at end of involved variable keys */
    std::vector<int>::const_iterator end() const
    {
        return keys_.end();
    }
    /**
    * @return the number of variables involved in this factor
    */
    int size() const
    {
        return keys_.size();
    }

    void push_back(int key)
    {
        keys_.push_back(key);
    }
    void clear()
    {
        keys_.clear();
    }

    Factor& operator=(const Factor& rObj)
    {
        this->keys_ = rObj.keys_;
        return *this;
    }

    /// @}


    /// @name Testable
    /// @{

protected:
    /// check equality
    bool equals(const Factor& other, double tol = 1e-9) const
    {
        return keys_ == other.keys_;
    };

    /// @}

public:
    /// @name Advanced Interface
    /// @{

    /** @return keys involved in this factor */
    std::vector<int>& keys()
    {
        return keys_;
    }

    /** Iterator at beginning of involved variable keys */
    std::vector<int>::iterator begin()
    {
        return keys_.begin();
    }

    /** Iterator at end of involved variable keys */
    std::vector<int>::iterator end()
    {
        return keys_.end();
    }

    void erase(std::vector<int>::iterator needtoerase)
    {
        keys_.erase(needtoerase);
    }
    void erase(int key)
    {
        //keys_.erase();
        keys_.erase(keys_.begin()+key);
    }
    bool empty()
    {
        return keys_.empty();
    }

    bool operator==(const Factor& other)const;
};

//}


#endif // FACTOR_H_INCLUDED
