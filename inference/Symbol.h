#ifndef SYMBOL_H
#define SYMBOL_H


/**
 * @file Symbol.h
 * @date
 * @author:
 */


#include <cstdint>
#include <string>

/**
 * Character and index key used in VectorValues, GaussianFactorGraph,
 * GaussianFactor, etc.  These keys are generated at runtime from TypedSymbol
 * keys when linearizing a nonlinear factor graph.  This key is not type
 * safe, so cannot be used with any Nonlinear* classes.
 */
class  Symbol
{
public:
    unsigned char c_;
    int j_;

public:

    /** Default constructor */
    Symbol() :
        c_(0), j_(0)
    {
    }

    /** Copy constructor */
    Symbol(const Symbol& key) :
        c_(key.c_), j_(key.j_)
    {
    }

    /** Constructor */
    Symbol(unsigned char c, int j) :
        c_(c), j_(j)
    {
    }

    /** Constructor that decodes an integer Key */
    Symbol(int key);

    /** return Key (integer) representation */
    int key() const;

    /** Cast to integer */
    int Key() const
    {
        return key();
    }

    /// Print
    //void print(const std::string& s = "") const;

    /// Check equality
    //bool equals(const Symbol& expected, double tol = 0.0) const;

    /** Retrieve key character */
    unsigned char chr() const
    {
        return c_;
    }

    /** Retrieve key index */
    int index() const
    {
        return j_;
    }

    /** Create a string from the key */
    std::string createstring() const
    {
        return std::string("%c%d",this->c_,this->j_);
    }

    /** Comparison for use in maps */
    bool operator<(const Symbol& comp) const
    {
        return c_ < comp.c_ || (comp.c_ == c_ && j_ < comp.j_);
    }

    /** Comparison for use in maps */
    bool operator==(const Symbol& comp) const
    {
        return comp.c_ == c_ && comp.j_ == j_;
    }

    /** Comparison for use in maps */
// bool operator==(Symbol comp) const {
//   return comp == (*this);
// }

    /** Comparison for use in maps */
    bool operator!=(const Symbol& comp) const
    {
        return comp.c_ != c_ || comp.j_ != j_;
    }

    /** Comparison for use in maps */
    bool operator!=(Symbol& comp) const
    {
        return comp != (*this);
    }

    /** Return a filter function that returns true when evaluated on a Key whose
     * character (when converted to a Symbol) matches \c c.  Use this with the
     * Values::filter() function to retrieve all key-value pairs with the
     * requested character.
     */
    //static boost::function<bool(Key)> ChrTest(unsigned char c);

//private:

    /** Serialization function */
    //friend class boost::serialization::access;
    //template<class ARCHIVE>
    //void serialize(ARCHIVE & ar, const unsigned int /*version*/) {
    //  ar & BOOST_SERIALIZATION_NVP(c_);
    //  ar & BOOST_SERIALIZATION_NVP(j_);
    //}
};

/** Create a symbol key from a character and index, i.e. x5. */
inline Symbol  symbol(unsigned char c, int j)
{
    return Symbol(c,j);
}

/** Return the character portion of a symbol key. */
inline unsigned char symbolChr(int key)
{
    return Symbol(key).chr();
}

/** Return the index portion of a symbol key. */
inline int symbolIndex(int key)
{
    return Symbol(key).index();
}

namespace symbol_shorthand
{
inline Symbol A(int j)
{
    return Symbol('a', j);
}
inline Symbol B(int j)
{
    return Symbol('b', j);
}
inline Symbol C(int j)
{
    return Symbol('c', j);
}
inline Symbol D(int j)
{
    return Symbol('d', j);
}
inline Symbol E(int j)
{
    return Symbol('e', j);
}
inline Symbol F(int j)
{
    return Symbol('f', j);
}
inline Symbol G(int j)
{
    return Symbol('g', j);
}
inline Symbol H(int j)
{
    return Symbol('h', j);
}
inline Symbol I(int j)
{
    return Symbol('i', j);
}
inline Symbol J(int j)
{
    return Symbol('j', j);
}
inline Symbol K(int j)
{
    return Symbol('k', j);
}
inline Symbol L(int j)
{
    return Symbol('l', j);
}
inline Symbol M(int j)
{
    return Symbol('m', j);
}
inline Symbol N(int j)
{
    return Symbol('n', j);
}
inline Symbol O(int j)
{
    return Symbol('o', j);
}
inline Symbol P(int j)
{
    return Symbol('p', j);
}
inline Symbol Q(int j)
{
    return Symbol('q', j);
}
inline Symbol R(int j)
{
    return Symbol('r', j);
}
inline Symbol S(int j)
{
    return Symbol('s', j);
}
inline Symbol T(int j)
{
    return Symbol('t', j);
}
inline Symbol U(int j)
{
    return Symbol('u', j);
}
inline Symbol V(int j)
{
    return Symbol('v', j);
}
inline Symbol W(int j)
{
    return Symbol('w', j);
}
inline Symbol X(int j)
{
    return Symbol('x', j);
}
inline Symbol Y(int j)
{
    return Symbol('y', j);
}
inline Symbol Z(int j)
{
    return Symbol('z', j);
}
}

#endif // SYMBOL_H
