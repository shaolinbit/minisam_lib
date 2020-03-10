#ifndef GAUSSIANLOCKMATRIX_H_INCLUDED
#define GAUSSIANLOCKMATRIX_H_INCLUDED

#include <stdio.h>
#include <array>
#include <stdexcept>
#include "Matrix.h"
#include "MatCal.h"
#include <iostream>



using namespace std;
namespace minisam
{

class GaussianBlockMatrix
{

public:
    bool BlockkindSB;                 /// if it is a square matrix
    int rowStart_;                    ///< Changes apparent matrix view, see main class comment.
    int rowEnd_;                      ///< Changes apparent matrix view, see main class comment.
    int blockStart_;                  ///< Changes apparent matrix view, see main class comment.
    minimatrix matrix_;

    vector<int> *variableColOffsets_; //Offsets Variable for Column

    ///constructors for initial the basic information of matrix
    GaussianBlockMatrix();
    GaussianBlockMatrix(bool blockkindsb);
    GaussianBlockMatrix(const vector<int> &dimensions, bool appendOneDimension, bool blockkindsb);
    GaussianBlockMatrix(vector<int>::iterator firstBlockDim, vector<int>::iterator lastBlockDim, bool appendOneDimension = false, bool blockkindsb = true);
    GaussianBlockMatrix(vector<int>::const_iterator firstBlockDim, vector<int>::const_iterator lastBlockDim, bool appendOneDimension = false, bool blockkindsb = true);
    GaussianBlockMatrix(const vector<int> &dimensions, const minimatrix& matrix,
                        bool appendOneDimension = false, bool blockkindsb = true);
    GaussianBlockMatrix(const GaussianBlockMatrix &rObj);
    GaussianBlockMatrix &operator=(const GaussianBlockMatrix &rObj);
    ~GaussianBlockMatrix();

    void setGaussianBlockMatrix(const GaussianBlockMatrix &rObj); ///initial this class by other class
    void assertInvariants() const;
    minimatrix& matrix();
    //// Get the basic information of the matrix
    int rows() const;    /// number of  Row
    int cols() const;    ////// number of Column
    int nBlocks() const; /// Block count

    /// Number of dimensions for variable on this diagonal block.
    int SgetDim(int block) const;

    minimatrix Sblock(int I, int J) const;


    /// Return the J'th diagonal block as a self adjoint view.
    //  minimatrix SdiagonalBlockNonconstant(int J);

    void SdiagonalBlock_dsyrk(const minimatrix& Ab_j,int J);

    minimatrix SdiagonalBlock(int J) const;

    /// Get the diagonal of the J'th diagonal block.
    minivector diagonal(int J) const;

    /// Get block above the diagonal (I, J).
    minimatrix SaboveDiagonalBlock(int I, int J) const;

    /// Return the square sub-matrix that contains blocks(i:j, i:j).
    minimatrix  selfadjointView(int I, int J) const;


    /// Return the square sub-matrix that contains blocks(i:j, i:j) as a
    /// triangular view.
    minimatrix triangularView(int I, int J) const;

    /// Get a range [i,j) from the matrix. Indices are in block units.
     minimatrix SconstaboveDiagonalRange(int i_startBlock, int i_endBlock, int j_startBlock, int j_endBlock) const;

    /// Get a range [i,j) from the matrix. Indices are in block units.
    minimatrix SaboveDiagonalRange(int i_startBlock, int i_endBlock, int j_startBlock, int j_endBlock);

    void SsetDiagonalBlock(int I, const minimatrix &xpr);
    void SsetOffDiagonalBlock(int I, int J, const minimatrix &xpr);
    /// Increment the diagonal block by the values in `xpr`. Only reads the upper
    /// triangular part of `xpr`.
    void SupdateDiagonalBlock(int I, const minimatrix& xpr);
    /// Update an off diagonal block.
    void SupdateOffDiagonalBlock(int I, int J, const minimatrix &xpr);
    /// @}
    /// @name Accessing the full matrix.
    /// @{

    /// Get self adjoint view.
    minimatrix selfadjointView();
    /// Get self adjoint view.
    const minimatrix selfadjointView() const;

    void SsetFullMatrix(const minimatrix&xpr);

    /// Set the upper active matrix zero.
    void UppersetZero();
    /// Set the entire active matrix zero.
    void setZero();


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
    GaussianBlockMatrix split(int nFrontals); //
    // V
    /** Construct from a container of the sizes of each vertical block. */
    GaussianBlockMatrix(vector<int> &dimensions, int height,
                        bool appendOneDimension, bool BlockkindSB);
    // V
    /** Construct from a container of the sizes of each vertical block and a
    * pre-prepared matrix. */
    GaussianBlockMatrix(vector<int> &dimensions, const minimatrix& matrix, bool appendOneDimension = false, bool BlockkindSB = false);

    /** Construct from iterator over the sizes of each vertical block. */
    GaussianBlockMatrix(vector<int>::iterator firstBlockDim, vector<int>::iterator lastBlockDim, int height, bool appendOneDimension = false, bool BlockkindSB = false);
    GaussianBlockMatrix(vector<int>::const_iterator firstBlockDim, vector<int>::const_iterator lastBlockDim, int height, bool appendOneDimension = false, bool BlockkindSB = false);

    /** Access a single block in the underlying matrix with read/write access
    * //operator()*/
    minimatrix VblockNonconst(int blocknum);
     minimatrix Vblock(int blocknum) const;
    /** access ranges of blocks at a time */
    minimatrix VrangeNonconst(int startBlock, int endBlock);

     minimatrix Vrange(int startBlock, int endBlock) const;

    void setVblock(int blocknum,const minimatrix& Ai);

    /** Return the full matrix, *not* including any portions excluded by
    * rowStart(), rowEnd(), and firstBlock() */
    minimatrix Vfull();
    /** Return the full matrix, *not* including any portions excluded by
    * rowStart(), rowEnd(), and firstBlock() */
    const minimatrix Vfull() const;

    void setVfull(const minimatrix& Ai);



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
     minimatrix matrix() const;

    /// Get an offset for a block index (in the active view).
    int Soffset(int block) const;

protected:
    int nOffsets() const;      /// Number of offsets in the full matrix.
    /// Get an arbitrary block from the matrix. Indices are in block units.
     minimatrix Sblock_(int iBlock, int jBlock, int blockRows = 1, int blockCols = 1) const;
    /// Get an arbitrary block from the matrix. Indices are in block units.
    minimatrix Sblock_nonconst(int iBlock, int jBlock, int blockRows = 1, int blockCols = 1);
    /// Get the full matrix as a block.
    const minimatrix Sfull() const;
    /// Get the full matrix as a block.
    minimatrix Sfullnonconst();

    /// Compute the indices into the underlying matrix for a given block.
    mini_int_vector calcIndices(int iBlock, int jBlock, int blockRows, int blockCols) const;

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

GaussianBlockMatrix LikeActiveViewOfVS(const GaussianBlockMatrix &other);
GaussianBlockMatrix LikeActiveViewOfSV(const GaussianBlockMatrix &rhs, int height);
} // namespace minisam

#endif // GAUSSIANLOCKMATRIX_H_INCLUDED

