#ifndef ISAM2DATA_H_INCLUDED
#define ISAM2DATA_H_INCLUDED
#include "../nonlinear/NonlinearFactorGraph.h"
namespace minisam
{
class ISAM2Data
{
public:
    /** The current linearization point */
    std::map<int,minimatrix*> theta_;
    std::map<int,minimatrix*> resulttheta_;
    std::map<int,minivector> delta_;

    mutable std::map<int,minivector> deltaNewton_; // Only used when using Dogleg - stores the Gauss-Newton update
    mutable std::map<int,minivector> RgProd_; // Only used when using Dogleg - stores R*g and is updated incrementally
    NonlinearFactorGraph nonlinearFactors_;
    /** The current linear factors, which are only updated as needed */
    mutable GaussianFactorGraph linearFactors_;
    VariableIndex variableIndex_;
    mutable std::set<int> deltaReplacedMask_;
    /** Set of variables that are involved with linear factors from marginalized
     * variables and thus cannot have their linearization points changed. */
    std::set<int> fixedVariables_;
    ISAM2Data();
    ~ISAM2Data();

    void clearfactors();
    void clearvalues();
};

};
#endif // ISAM2DATA_H_INCLUDED
