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
    bool BlockkindSB;
    Eigen::MatrixXd* matrix_; ///< The full matrix
    std::vector<int>*
    variableColOffsets_; ///< the starting columns of each block (0-based)

    int rowStart_;   ///< Changes apparent matrix view, see main class comment.
    int rowEnd_;     ///< Changes apparent matrix view, see main class comment.
    int blockStart_; ///< Changes apparent matrix view, see main class comment.
public:
    /** Construct an empty VerticalBlockMatrix */
    SVBlockMatrix();

    ~SVBlockMatrix();
    SVBlockMatrix(bool blockkindsb);

    /** Non-const access to full matrix (*including* any portions excluded by
    * rowStart(), rowEnd(), and firstBlock()) */
    Eigen::MatrixXd &matrix();

    // S
    SVBlockMatrix(const std::vector<int> &dimensions,
                  bool appendOneDimension, bool blockkindsb);
    // S

    /// Construct from iterator over the sizes of each vertical block.
    SVBlockMatrix(std::vector<int>::iterator firstBlockDim,
                  std::vector<int>::iterator lastBlockDim,
                  bool appendOneDimension = false, bool blockkindsb = true);

    SVBlockMatrix(std::vector<int>::const_iterator firstBlockDim,
                  std::vector<int>::const_iterator lastBlockDim,
                  bool appendOneDimension = false, bool blockkindsb = true);

    // S
    /// Construct from a container of the sizes of each vertical block and a
    /// pre-prepared matrix.
    SVBlockMatrix(const std::vector<int> &dimensions,
                  const Eigen::MatrixXd &matrix, bool appendOneDimension = false,
                  bool blockkindsb = true);

    SVBlockMatrix &operator=(const SVBlockMatrix &rObj);//

    SVBlockMatrix(const SVBlockMatrix &rObj);
    void setSVBlockMatrix(const SVBlockMatrix &rObj);

    void assertInvariants() const;

    /// Row size
    int rows() const;

    /// Column size
    int cols() const;

    /// Block count
    int nBlocks() const;

    // S
    /// Number of dimensions for variable on this diagonal block.
    int SgetDim(int block) const;

    // S
    /// @name Block getter methods.
    /// @{

    /// Get a copy of a block (anywhere in the matrix).
    /// This method makes a copy - use the methods below if performance is
    /// critical.
    Eigen::MatrixXd Sblock(int I, int J) const;

    // S
    /// Return the J'th diagonal block as a self adjoint view.
    Eigen::SelfAdjointView<Eigen::Block<Eigen::MatrixXd>, Eigen::Upper>
    SdiagonalBlock(int J);

    // S
    /// Return the J'th diagonal block as a self adjoint view.
    Eigen::SelfAdjointView<Eigen::Block<const Eigen::MatrixXd>, Eigen::Upper>
    SdiagonalBlock(int J) const;

    // S
    /// Get the diagonal of the J'th diagonal block.
    Eigen::VectorXd diagonal(int J) const;
    // S
    /// Get block above the diagonal (I, J).
    Eigen::Block<const Eigen::MatrixXd> SaboveDiagonalBlock(int I, int J) const;

    // S
    /// Return the square sub-matrix that contains blocks(i:j, i:j).
    Eigen::SelfAdjointView<Eigen::Block<const Eigen::MatrixXd>, Eigen::Upper>
    selfadjointView(int I, int J) const;

    // S
    /// Return the square sub-matrix that contains blocks(i:j, i:j) as a
    /// triangular view.
    Eigen::TriangularView<Eigen::Block<const Eigen::MatrixXd>, Eigen::Upper>
    triangularView(int I, int J) const;

    // S
    /// Get a range [i,j) from the matrix. Indices are in block units.
    Eigen::Block<const Eigen::MatrixXd> SconstaboveDiagonalRange(
        int i_startBlock, int i_endBlock, int j_startBlock,
        int j_endBlock) const;

    Eigen::Block<const Eigen::MatrixXd> SaboveDiagonalRange(
        int i_startBlock, int i_endBlock, int j_startBlock,
        int j_endBlock) const;

    // S
    /// Get a range [i,j) from the matrix. Indices are in block units.
    Eigen::Block<Eigen::MatrixXd> SaboveDiagonalRange(int i_startBlock,
            int i_endBlock,
            int j_startBlock,
            int j_endBlock);

    /// @}
    /// @name Block setter methods.
    /// @{

    // S
    /// Set a diagonal block. Only the upper triangular portion of `xpr` is
    /// evaluated.
    void SsetDiagonalBlock(int I, const Eigen::MatrixXd &xpr);
    // S
    /// Set an off-diagonal block. Only the upper triangular portion of `xpr` is
    /// evaluated.
    void SsetOffDiagonalBlock(int I, int J, const Eigen::MatrixXd &xpr);
    // S
    /// Increment the diagonal block by the values in `xpr`. Only reads the upper
    /// triangular part of `xpr`.
    void SupdateDiagonalBlock(int I, const Eigen::MatrixXd &xpr);
    // S

    /// Update an off diagonal block.
    /// NOTE(emmett): This assumes noalias().
    void SupdateOffDiagonalBlock(int I, int J, const Eigen::MatrixXd &xpr);
    /// @}
    /// @name Accessing the full matrix.
    /// @{

    // S
    /// Get self adjoint view.
    Eigen::SelfAdjointView<Eigen::Block<Eigen::MatrixXd>, Eigen::Upper>
    selfadjointView();

    // S
    /// Get self adjoint view.
    Eigen::SelfAdjointView<Eigen::Block<const Eigen::MatrixXd>, Eigen::Upper>
    selfadjointView() const;

    // S
    /// Set the entire active matrix. Only reads the upper triangular part of
    /// `xpr`.
    void SsetFullMatrix(const Eigen::MatrixXd &xpr);

    // S
    /// Set the entire active matrix zero.
    void SsetZero();


    /// @}

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

    // StoV
    SVBlockMatrix split(int nFrontals); //

public: // V
    // V
    /** Construct from a container of the sizes of each vertical block. */
    SVBlockMatrix(std::vector<int> &dimensions, int height,
                  bool appendOneDimension, bool BlockkindSB);
    // V
    /** Construct from a container of the sizes of each vertical block and a
    * pre-prepared matrix. */
    SVBlockMatrix(std::vector<int> &dimensions,
                  const Eigen::MatrixBase<Eigen::MatrixXd> &matrix,
                  bool appendOneDimension = false, bool BlockkindSB = false);

    /** Construct from iterator over the sizes of each vertical block. */
    SVBlockMatrix(std::vector<int>::iterator firstBlockDim,
                  std::vector<int>::iterator lastBlockDim, int height,
                  bool appendOneDimension = false, bool BlockkindSB = false);
    SVBlockMatrix(std::vector<int>::const_iterator firstBlockDim,
                  std::vector<int>::const_iterator lastBlockDim, int height,
                  bool appendOneDimension = false, bool BlockkindSB = false);

    /** Access a single block in the underlying matrix with read/write access
    * //operator()*/
    Eigen::Block<Eigen::MatrixXd> Vblock(int blocknum);

    Eigen::Block<const  Eigen::MatrixXd> Vblock(int blocknum)  const;

    /** access ranges of blocks at a time */
    Eigen::Block<Eigen::MatrixXd> Vrange(int startBlock, int endBlock);
    Eigen::Block<const Eigen::MatrixXd> Vrange(int startBlock,
            int endBlock) const;
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


public:
    // S
    /// Get an offset for a block index (in the active view).
    int Soffset(int block) const;

protected:
    // S
    /// Number of offsets in the full matrix.
    int nOffsets() const;

    // S
    /// Number of actual blocks in the full matrix.
    int nActualBlocks() const;

    // S
    /// Get an arbitrary block from the matrix. Indices are in block units.
    Eigen::Block<const Eigen::MatrixXd> Sblock_(int iBlock, int jBlock,
            int blockRows = 1,
            int blockCols = 1) const;

    // S
    /// Get an arbitrary block from the matrix. Indices are in block units.
    Eigen::Block<Eigen::MatrixXd> Sblock_(int iBlock, int jBlock,
                                          int blockRows = 1, int blockCols = 1);

    // S
    /// Get the full matrix as a block.
    Eigen::Block<const Eigen::MatrixXd> Sfull() const;

    // S
    /// Get the full matrix as a block.
    Eigen::Block<Eigen::MatrixXd> Sfull();

    /// Compute the indices into the underlying matrix for a given block.
    std::array<int, 4> calcIndices(int iBlock, int jBlock, int blockRows,
                                   int blockCols) const;
    void SfillOffsets(std::vector<int>::const_iterator firstBlockDim,
                      std::vector<int>::const_iterator lastBlockDim,
                      bool appendOneDimension);

    /// V
    void checkBlock(int block) const;

    void VfillOffsets(std::vector<int>::iterator firstBlockDim,
                      std::vector<int>::iterator lastBlockDim,
                      bool appendOneDimension);
    void VfillOffsets(std::vector<int>::const_iterator firstBlockDim,
                      std::vector<int>::const_iterator lastBlockDim,
                      bool appendOneDimension);
};

SVBlockMatrix  LikeActiveViewOfVS(const SVBlockMatrix& other);
SVBlockMatrix LikeActiveViewOfSV(const SVBlockMatrix& rhs, int height);
}

#endif // SVBLOCKMATRIX_H_INCLUDED
