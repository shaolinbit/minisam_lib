#ifndef SVBLOCKMATRIX_H_INCLUDED
#define SVBLOCKMATRIX_H_INCLUDED

#include <stdio.h>
#include <array>
#include <stdexcept>
#include "../base/Matrix.h"
#include "../base/MatCal.h"

#include <iostream>
/*

Combine two kinds of blockmatrix, VerticalBlockMatrix and SymmetricBlockMatrix, to one.

*/
/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    VerticalBlockMatrix.h
 * @brief   A matrix with column blocks of pre-defined sizes.  Used in JacobianFactor and
 *          GaussianConditional.
 * @author  Richard Roberts
 * @date    Sep 18, 2010 */

/**
* @file    SymmetricBlockMatrix.h
* @brief   Access to matrices via blocks of pre-defined sizes.  Used in GaussianFactor and GaussianConditional.
* @author  Richard Roberts
* @date    Sep 18, 2010
*/



using namespace std;
namespace minisam
{

class SVBlockMatrix
{

public:
    //Date member,private member
    bool BlockkindSB;                 /// if it is a square matrix
    int rowStart_;                    ///< Changes apparent matrix view, see main class comment.
    int rowEnd_;                      ///< Changes apparent matrix view, see main class comment.
    int blockStart_;                  ///< Changes apparent matrix view, see main class comment.
    Eigen::MatrixXd *matrix_;         ///

    vector<int> *variableColOffsets_; //Offsets Variable for Column

    ///constructors for initial the basic information of matrix
    SVBlockMatrix();                                                                         //checked already                                                             //default constructors                                                                   ///构造函数
    SVBlockMatrix(bool blockkindsb);                                                         //checked already                                                     ///构造函数
    SVBlockMatrix(const vector<int> &dimensions, bool appendOneDimension, bool blockkindsb); ///重载构造函数 //already                                                                                   //复制构造函数
    SVBlockMatrix(vector<int>::iterator firstBlockDim, vector<int>::iterator lastBlockDim, bool appendOneDimension = false, bool blockkindsb = true);
    SVBlockMatrix(vector<int>::const_iterator firstBlockDim, vector<int>::const_iterator lastBlockDim, bool appendOneDimension = false, bool blockkindsb = true);
    SVBlockMatrix(const vector<int> &dimensions, const Eigen::MatrixXd &matrix, bool appendOneDimension = false, bool blockkindsb = true);
    SVBlockMatrix(const SVBlockMatrix &rObj);
    SVBlockMatrix &operator=(const SVBlockMatrix &rObj); ///赋值符号重载
    ~SVBlockMatrix();                                    ///析构函数

    void setSVBlockMatrix(const SVBlockMatrix &rObj); ///initial this class by other class
    void assertInvariants() const;                    //this function is for error checking
    Eigen::MatrixXd &matrix();
    //// Get the basic information of the matrix
    int rows() const;    /// number of  Row
    int cols() const;    ////// number of Column
    int nBlocks() const; /// Block count

    /// Number of dimensions for variable on this diagonal block.
    int SgetDim(int block) const;
    Eigen::MatrixXd Sblock(int I, int J) const;
    /// Return the J'th diagonal block as a self adjoint view.
    Eigen::SelfAdjointView<Eigen::Block<Eigen::MatrixXd>, Eigen::Upper> SdiagonalBlock(int J);
    Eigen::SelfAdjointView<Eigen::Block<const Eigen::MatrixXd>, Eigen::Upper> SdiagonalBlock(int J) const;

    /// Get the diagonal of the J'th diagonal block.
    Eigen::VectorXd diagonal(int J) const;
    /// Get block above the diagonal (I, J).
    Eigen::Block<const Eigen::MatrixXd> SaboveDiagonalBlock(int I, int J) const;

    /// Return the square sub-matrix that contains blocks(i:j, i:j).
    Eigen::SelfAdjointView<Eigen::Block<const Eigen::MatrixXd>, Eigen::Upper>
    selfadjointView(int I, int J) const;

    /// Return the square sub-matrix that contains blocks(i:j, i:j) as a
    /// triangular view.
    Eigen::TriangularView<Eigen::Block<const Eigen::MatrixXd>, Eigen::Upper>
    triangularView(int I, int J) const;

    /// Get a range [i,j) from the matrix. Indices are in block units.
    Eigen::Block<const Eigen::MatrixXd> SconstaboveDiagonalRange(int i_startBlock, int i_endBlock, int j_startBlock, int j_endBlock) const;
    Eigen::Block<const Eigen::MatrixXd> SaboveDiagonalRange(int i_startBlock, int i_endBlock, int j_startBlock, int j_endBlock) const;

    /// Get a range [i,j) from the matrix. Indices are in block units.
    Eigen::Block<Eigen::MatrixXd> SaboveDiagonalRange(int i_startBlock, int i_endBlock, int j_startBlock, int j_endBlock);

    void SsetDiagonalBlock(int I, const Eigen::MatrixXd &xpr);
    void SsetOffDiagonalBlock(int I, int J, const Eigen::MatrixXd &xpr);
    /// Increment the diagonal block by the values in `xpr`. Only reads the upper
    /// triangular part of `xpr`.
    void SupdateDiagonalBlock(int I, const Eigen::MatrixXd &xpr);

    /// Update an off diagonal block.
    /// NOTE(emmett): This assumes noalias().
    void SupdateOffDiagonalBlock(int I, int J, const Eigen::MatrixXd &xpr);
    /// @}
    /// @name Accessing the full matrix.
    /// @{

    /// Get self adjoint view.
    Eigen::SelfAdjointView<Eigen::Block<Eigen::MatrixXd>, Eigen::Upper> selfadjointView();
    /// Get self adjoint view.
    Eigen::SelfAdjointView<Eigen::Block<const Eigen::MatrixXd>, Eigen::Upper> selfadjointView() const;

    void SsetFullMatrix(const Eigen::MatrixXd &xpr);

    /// Set the entire active matrix zero.
    void SsetZero();

    /// Retrieve or modify the first logical block, i.e. the block referenced by
    /// block index 0. Blocks before it will be inaccessible, except by accessing
    /// the underlying matrix using matrix().
    int &blockStart();

    /// Retrieve the first logical block, i.e. the block referenced by block index
    /// 0. Blocks before it will be inaccessible, except by accessing the
    /// underlying matrix using matrix().
    int blockStart() const;

    /**
    * After partial Cholesky, we can optionally split off R and Sd, to be
    * interpreted as a GaussianConditional |R*x1 + S*x2 - d]^2. We leave the
    * symmetric lower block L in place, and adjust block_start so now *this
    * refers to it.
    */
    SVBlockMatrix
    split(int nFrontals); //
    // V
    /** Construct from a container of the sizes of each vertical block. */
    SVBlockMatrix(vector<int> &dimensions, int height,
                  bool appendOneDimension, bool BlockkindSB);
    // V
    /** Construct from a container of the sizes of each vertical block and a
    * pre-prepared matrix. */
    SVBlockMatrix(vector<int> &dimensions, const Eigen::MatrixBase<Eigen::MatrixXd> &matrix, bool appendOneDimension = false, bool BlockkindSB = false);

    /** Construct from iterator over the sizes of each vertical block. */
    SVBlockMatrix(vector<int>::iterator firstBlockDim, vector<int>::iterator lastBlockDim, int height, bool appendOneDimension = false, bool BlockkindSB = false);
    SVBlockMatrix(vector<int>::const_iterator firstBlockDim, vector<int>::const_iterator lastBlockDim, int height, bool appendOneDimension = false, bool BlockkindSB = false);

    /** Access a single block in the underlying matrix with read/write access
    * //operator()*/
    Eigen::Block<Eigen::MatrixXd> Vblock(int blocknum);
    Eigen::Block<const Eigen::MatrixXd> Vblock(int blocknum) const;

    /** access ranges of blocks at a time */
    Eigen::Block<Eigen::MatrixXd> Vrange(int startBlock, int endBlock);
    Eigen::Block<const Eigen::MatrixXd> Vrange(int startBlock, int endBlock) const;
    /** Return the full matrix, *not* including any portions excluded by
    * rowStart(), rowEnd(), and firstBlock() */
    Eigen::Block<Eigen::MatrixXd> Vfull();

    /** Return the full matrix, *not* including any portions excluded by
    * rowStart(), rowEnd(), and firstBlock() */
    const Eigen::Block<const Eigen::MatrixXd> Vfull() const;

    int Voffset(int block) const;
    /** Get the apparent first row of the underlying matrix for all operations */
    const int &rowStart() const;

    /** Get or set the apparent first row of the underlying matrix for all
    * operations */
    int &rowStart();

    /** Get the apparent last row (exclusive, i.e. rows() == rowEnd() -
    * rowStart()) of the underlying matrix for all operations */
    const int &rowEnd() const;

    /** Get or set the apparent last row (exclusive, i.e. rows() == rowEnd() -
    * rowStart()) of the underlying matrix for all operations */
    int &rowEnd();

    /** Get the apparent first block for all operations */
    const int &firstBlock() const;

    /** Get or set the apparent first block for all operations */
    int &firstBlock();

    /** Access to full matrix (*including* any portions excluded by rowStart(),
    * rowEnd(), and firstBlock()) */
    const Eigen::MatrixXd &matrix() const;
    /// Get an offset for a block index (in the active view).
    int Soffset(int block) const;

protected:
    int nOffsets() const;      /// Number of offsets in the full matrix.
    int nActualBlocks() const; /// Number of actual blocks in the full matrix.
    /// Get an arbitrary block from the matrix. Indices are in block units.
    Eigen::Block<const Eigen::MatrixXd> Sblock_(int iBlock, int jBlock, int blockRows = 1, int blockCols = 1) const;

    /// Get an arbitrary block from the matrix. Indices are in block units.
    Eigen::Block<Eigen::MatrixXd> Sblock_(int iBlock, int jBlock, int blockRows = 1, int blockCols = 1);
    /// Get the full matrix as a block.
    Eigen::Block<const Eigen::MatrixXd> Sfull() const;
    /// Get the full matrix as a block.
    Eigen::Block<Eigen::MatrixXd> Sfull();
    /// Compute the indices into the underlying matrix for a given block.
    array<int, 4> calcIndices(int iBlock, int jBlock, int blockRows, int blockCols) const;
    void SfillOffsets(vector<int>::const_iterator firstBlockDim, vector<int>::const_iterator lastBlockDim, bool appendOneDimension);

    /// V
    void checkBlock(int block) const;

    void VfillOffsets(vector<int>::iterator firstBlockDim,
                      vector<int>::iterator lastBlockDim,
                      bool appendOneDimension);
    void VfillOffsets(vector<int>::const_iterator firstBlockDim,
                      vector<int>::const_iterator lastBlockDim,
                      bool appendOneDimension);
};

SVBlockMatrix LikeActiveViewOfVS(const SVBlockMatrix &other);           // 全局变量，类的实列化
SVBlockMatrix LikeActiveViewOfSV(const SVBlockMatrix &rhs, int height); //全局变量，全局变量类的实际例化
} // namespace minisam

#endif // SVBLOCKMATRIX_H_INCLUDED

