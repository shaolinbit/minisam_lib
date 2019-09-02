#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    Matrix.h
 * @brief   typedef and functions to augment Eigen's MatrixXd
 * @author  Christian Potthast
 * @author  Kai Ni
 * @author  Frank Dellaert
 * @author  Alex Cunningham
 * @author  Alex Hagiopol
 */


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



namespace minisam
{
///////////////////////
void zeroBelowDiagonal(Eigen::MatrixXd& A, size_t cols=0);
double norm2d(const Eigen::Vector3d& p);
double norm2d(const Eigen::Vector2d& p,  Eigen::MatrixXd* H);
void vector_scale_inplace(const Eigen::VectorXd& v, Eigen::MatrixXd& A, bool inf_mask = false); // row
Eigen::MatrixXd vector_scale(const Eigen::VectorXd& v, const Eigen::MatrixXd& A, bool inf_mask = false); // row
Eigen::MatrixXd vector_scale(const Eigen::MatrixXd& A, const Eigen::VectorXd& v, bool inf_mask = false); // column
double distance2(const Eigen::Vector2d& p, const Eigen::Vector2d& q);


inline Eigen::Matrix3d skewSymmetric(double wx, double wy, double wz)
{
    Eigen::Matrix3d m3d;
    m3d<<0.0, -wz, +wy, +wz, 0.0, -wx, -wy, +wx, 0.0;
    return m3d;
}

inline Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d& w)
{
    return skewSymmetric(w(0), w(1), w(2));
}
bool equal_with_abs_tol(const Eigen::MatrixXd A, const Eigen::MatrixXd B, double tol=1e-9);
bool equal_with_abs_tol(const Eigen::Matrix3d A, const Eigen::Matrix3d B, double tol=1e-9);

Eigen::VectorXd backSubstituteUpper(const Eigen::MatrixXd& U, const Eigen::VectorXd& b, bool unit=false);

void inplace_QR(Eigen::MatrixXd& A);

bool assert_equal(const Eigen::VectorXd& expected, const Eigen::VectorXd& actual, double tol);
bool equal_with_abs_tol(const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2, double tol);
bool equal_with_abs_tol(const Eigen::Vector3d& vec1, const Eigen::Vector3d& vec2, double tol=1e-9);


std::pair<int,bool> choleskyCareful(Eigen::MatrixXd& ATA, int order = -1);

bool choleskyPartial(Eigen::MatrixXd* ABC, int nFrontal, int topleft=0);

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
};
#endif // MATRIX_H_INCLUDED
