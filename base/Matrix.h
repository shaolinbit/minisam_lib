#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <Eigen/Geometry>
#include <vector>
#include <map>
#include <list>
#include "../gmfconfig.h"


#ifndef M_PI
        #define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
        #define M_PI_2 M_PI*0.5
#endif

std::vector<int> SplitVectorIntFrontal(std::vector<int> keys,int frontals);


std::vector<int> SplitVectorIntParents(std::vector<int> keys,int frontals);

double dot3d(const Eigen::Vector3d &p, const Eigen::Vector3d &q,
             Eigen::MatrixXd* H1,
             Eigen::MatrixXd* H2);
double dot3d(const Eigen::Vector3d &p, const Eigen::Vector3d &q);

double norm3d(const Eigen::Vector3d &p);
double norm3d(const Eigen::Vector3d &p, Eigen::MatrixXd* H);
double norm2d(const Eigen::Vector3d& p);
double norm2d(const Eigen::Vector2d& p,  Eigen::MatrixXd* H);

double distance2(const Eigen::Vector2d& p, const Eigen::Vector2d& q);

Eigen::Vector3d cross3d(const Eigen::Vector3d &p,
                        const Eigen::Vector3d &q,
                        Eigen::Matrix3d* H1,
                        Eigen::Matrix3d* H2);
Eigen::Vector3d cross3d(const Eigen::Vector3d &p,
                        const Eigen::Vector3d &q,
                        Eigen::Matrix3d* H1);
Eigen::Vector3d cross3d(const Eigen::Vector3d &p,
                        const Eigen::Vector3d &q);
Eigen::Vector3d normalize3d(const Eigen::Vector3d &p,
                            Eigen::Matrix3d* H);
Eigen::Vector3d normalize3d(const Eigen::Vector3d &p);
bool equal_with_abs_tol(const Eigen::MatrixXd A, const Eigen::MatrixXd B, double tol=1e-9);
bool equal_with_abs_tol(const Eigen::Matrix3d A, const Eigen::Matrix3d B, double tol=1e-9);

Eigen::VectorXd backSubstituteUpper(const Eigen::MatrixXd& U, const Eigen::VectorXd& b, bool unit=false);

void inplace_QR(Eigen::MatrixXd& A);

void zeroBelowDiagonal(Eigen::MatrixXd& A, size_t cols=0);

void vector_scale_inplace(const Eigen::VectorXd& v, Eigen::MatrixXd& A, bool inf_mask = false); // row
Eigen::MatrixXd vector_scale(const Eigen::VectorXd& v, const Eigen::MatrixXd& A, bool inf_mask = false); // row
Eigen::MatrixXd vector_scale(const Eigen::MatrixXd& A, const Eigen::VectorXd& v, bool inf_mask = false); // column


/**
 * "Careful" Cholesky computes the positive square-root of a positive symmetric
 * semi-definite matrix (i.e. that may be rank-deficient).  Unlike standard
 * Cholesky, the square-root factor may have all-zero rows for free variables.
 *
 * Additionally, this function returns the index of the row after the last
 * non-zero row in the computed factor, so that it may be truncated to an
 * upper-trapazoidal matrix.
 *
 * The second element of the return value is \c true if the matrix was factored
 * successfully, or \c false if it was non-positive-semidefinite (i.e.
 * indefinite or negative-(semi-)definite.
 *
 * Note that this returned index is the rank of the matrix if and only if all
 * of the zero-rows of the factor occur after any non-zero rows.  This is
 * (always?) the case during elimination of a fully-constrained least-squares
 * problem.
 *
 * The optional order argument specifies the size of the square upper-left
 * submatrix to operate on, ignoring the rest of the matrix.
 *
 *
 */
std::pair<int,bool> choleskyCareful(Eigen::MatrixXd& ATA, int order = -1);

/**
 * Partial Cholesky computes a factor [R S  such that [R' 0  [R S  = [A  B
 *                                     0 L]            S' I]  0 L]    B' C].
 * The input to this function is the matrix ABC = [A  B], and the parameter
 *                                                [B' C]
 * nFrontal determines the split between A, B, and C, with A being of size
 * nFrontal x nFrontal.
 *
 * if non-zero, factorization proceeds in bottom-right corner starting at topleft
 *
 * @return \c true if the decomposition is successful, \c false if \c A was
 * not positive-definite.
 */
bool choleskyPartial(Eigen::MatrixXd* ABC, int nFrontal, int topleft=0);

Eigen::VectorXd VectorValuesgetvector(const std::map<int,Eigen::VectorXd>& VectorValues,const std::vector<int>& keys);


template <class TPFactor>
std::vector<int>& getvectorint_fromiterator(TPFactor& nrParents_);


double VectorValuesDot(std::map<int,Eigen::VectorXd>& vectorvalues);

double VectorValuesDot(const std::map<int,Eigen::VectorXd>& vectorvalues1,const std::map<int,Eigen::VectorXd>& vectorvalues2);

//template<class TPFactor>
void VectorValuesScaleMultiply(std::map<int,Eigen::VectorXd> *vectorvalues,double a);

void VectorValuesaxpy(std::map<int,Eigen::VectorXd> *y,
                      const std::map<int,Eigen::VectorXd>& x,double a);

double ListVectorDot(std::list<Eigen::VectorXd>& Errors);


std::map<int,Eigen::VectorXd> VectorValuesZero(std::map<int,Eigen::VectorXd> vvz);





//void updateHessian(const std::vector<int>& keys, SymmetricBlockMatrix* info) ;
std::list<Eigen::VectorXd> nullstlevx();

std::map<int,Eigen::VectorXd> VectorValuesRetract(const std::map<int,Eigen::VectorXd>& vectorvalues1,
        const std::map<int,Eigen::VectorXd>& vectorvalues2);
std::map<int,Eigen::VectorXd> DVectorValuesRetract(const std::map<int,Eigen::VectorXd>& vectorvalues1,
        const std::map<int,Eigen::VectorXd>& vectorvalues2);

std::map<int,Eigen::VectorXd> VectorValuesSub(const std::map<int,Eigen::VectorXd>& vectorvalues1,
        const std::map<int,Eigen::VectorXd>& vectorvalues2);

std::map<int,Eigen::VectorXd> VectorValuesAdd(const std::map<int,Eigen::VectorXd>& vectorvalues1,
        const std::map<int,Eigen::VectorXd>& vectorvalues2);



void VectorValuesUpdate(std::map<int,Eigen::VectorXd>* vectorvalues1,const std::map<int,Eigen::VectorXd>& vectorvalues2);

double VectorValuesNorm(std::map<int,Eigen::VectorXd> vvz);


double VectorValuesSquareNorm(std::map<int,Eigen::VectorXd> vvz);

inline Eigen::Matrix3d skewSymmetric(double wx, double wy, double wz)
{
    Eigen::Matrix3d m3d;
    m3d<<0.0, -wz, +wy, +wz, 0.0, -wx, -wy, +wx, 0.0;
    return m3d;
    //return (Matrix3() << 0.0, -wz, +wy, +wz, 0.0, -wx, -wy, +wx, 0.0).finished();
}

//template <class Derived>
inline Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d w)
{
    return skewSymmetric(w(0), w(1), w(2));
}
bool assert_equal(const Eigen::VectorXd& expected, const Eigen::VectorXd& actual, double tol);
bool equal_with_abs_tol(const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2, double tol);
bool equal_with_abs_tol(const Eigen::Vector3d& vec1, const Eigen::Vector3d& vec2, double tol=1e-9);

std::vector<int> getvectorfrom3keys(int key1,int key2,int key3);
std::vector<int> getvectorfrom5keys(int key1,int key2,int key3,int key4,int key5);
std::vector<int> getvectorfrom6keys(int key1,int key2,int key3,int key4,int key5,int key6);

Eigen::VectorXd ComposeTwoVec3(const Eigen::Vector3d& v1,const Eigen::Vector3d& v2);

Eigen::VectorXd Vec3ToVecX3(const Eigen::Vector3d& v);

Eigen::MatrixXd Mat3ToMatX3(const Eigen::Matrix3d& v);

struct GraphvizFormatting
{
    enum Axis { X, Y, Z, NEGX, NEGY, NEGZ }; ///< World axes to be assigned to paper axes
    Axis paperHorizontalAxis; ///< The world axis assigned to the horizontal paper axis
    Axis paperVerticalAxis; ///< The world axis assigned to the vertical paper axis
    double figureWidthInches; ///< The figure width on paper in inches
    double figureHeightInches; ///< The figure height on paper in inches
    double scale; ///< Scale all positions to reduce / increase density
    bool mergeSimilarFactors; ///< Merge multiple factors that have the same connectivity
    bool plotFactorPoints; ///< Plots each factor as a dot between the variables
    bool connectKeysToFactor; ///< Draw a line from each key within a factor to the dot of the factor
    bool binaryEdges; ///< just use non-dotted edges for binary factors
    std::map<int, Eigen::Vector2d> factorPositions; ///< (optional for each factor) Manually specify factor "dot" positions.
    /// Default constructor sets up robot coordinates.  Paper horizontal is robot Y,
    /// paper vertical is robot X.  Default figure size of 5x5 in.
    GraphvizFormatting() :
        paperHorizontalAxis(Y), paperVerticalAxis(X),
        figureWidthInches(5), figureHeightInches(5), scale(1),
        mergeSimilarFactors(false), plotFactorPoints(true),
        connectKeysToFactor(true), binaryEdges(true) {}
};

#endif // MATRIX_H_INCLUDED
