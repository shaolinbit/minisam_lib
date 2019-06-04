#ifndef SVBLOCKMATRIX_H_INCLUDED
#define SVBLOCKMATRIX_H_INCLUDED

#include <stdio.h>
#include <array>
#include <stdexcept>
#include "../base/Matrix.h"

////debug include
#include <iostream>
using namespace std;
////debug include

class SVBlockMatrix
{
public:
    //   typedef VerticalBlockMatrix This;
    //   typedef Eigen::Block<Matrix> Block;
    //  typedef Eigen::Block<const Matrix> constBlock;

    bool BlockkindSB;
    Eigen::MatrixXd* matrix_; ///< The full matrix
    std::vector<int>*
    variableColOffsets_; ///< the starting columns of each block (0-based)

    int rowStart_;   ///< Changes apparent matrix view, see main class comment.
    int rowEnd_;     ///< Changes apparent matrix view, see main class comment.
    int blockStart_; ///< Changes apparent matrix view, see main class comment.
public:
    /** Construct an empty VerticalBlockMatrix */
    SVBlockMatrix()
        : rowStart_(0), rowEnd_(0), blockStart_(0), BlockkindSB(true)
    {
       // variableColOffsets_.push_back(0);
       variableColOffsets_=new std::vector<int>();
       variableColOffsets_->push_back(0);
       matrix_=new Eigen::MatrixXd();
        assertInvariants();
    }

     ~SVBlockMatrix()
    {
      delete variableColOffsets_;
      variableColOffsets_=NULL;
      delete matrix_;
      matrix_=NULL;
    }

    SVBlockMatrix(bool blockkindsb)
        : rowStart_(0), rowEnd_(0), blockStart_(0), BlockkindSB(blockkindsb)
    {

        variableColOffsets_=new std::vector<int>();
        variableColOffsets_->push_back(0);
         matrix_=new Eigen::MatrixXd();
        assertInvariants();
    }

    /** Non-const access to full matrix (*including* any portions excluded by
    * rowStart(), rowEnd(), and firstBlock()) */
    Eigen::MatrixXd &matrix()
    {
        return *matrix_;
    }

    // S
    SVBlockMatrix(const std::vector<int> &dimensions,
                  bool appendOneDimension, bool blockkindsb)
        : blockStart_(0), rowEnd_(0), rowStart_(0), BlockkindSB(blockkindsb)
    {
        variableColOffsets_=new std::vector<int>();
        matrix_=new Eigen::MatrixXd();
        SfillOffsets(dimensions.begin(), dimensions.end(), appendOneDimension);
        matrix_->resize(variableColOffsets_->back(), variableColOffsets_->back());
        assertInvariants();
    }
    // S

    /// Construct from iterator over the sizes of each vertical block.
    // template<typename ITERATOR>
    SVBlockMatrix(std::vector<int>::iterator firstBlockDim,
                  std::vector<int>::iterator lastBlockDim,
                  bool appendOneDimension = false, bool blockkindsb = true)
        : blockStart_(0), rowEnd_(0), rowStart_(0), BlockkindSB(blockkindsb)
    {
        variableColOffsets_=new std::vector<int>();
        matrix_=new Eigen::MatrixXd();
        SfillOffsets(firstBlockDim, lastBlockDim, appendOneDimension);
        matrix_->resize(variableColOffsets_->back(), variableColOffsets_->back());
        assertInvariants();
    }

    SVBlockMatrix(std::vector<int>::const_iterator firstBlockDim,
                  std::vector<int>::const_iterator lastBlockDim,
                  bool appendOneDimension = false, bool blockkindsb = true)
        : blockStart_(0), rowEnd_(0), rowStart_(0), BlockkindSB(blockkindsb)
    {
        variableColOffsets_=new std::vector<int>();
        matrix_=new Eigen::MatrixXd();
        SfillOffsets(firstBlockDim, lastBlockDim, appendOneDimension);
        matrix_->resize(variableColOffsets_->back(), variableColOffsets_->back());
        assertInvariants();
    }

    // S
    /// Construct from a container of the sizes of each vertical block and a
    /// pre-prepared matrix.
    // template<typename CONTAINER>
    SVBlockMatrix(const std::vector<int> &dimensions,
                  const Eigen::MatrixXd &matrix, bool appendOneDimension = false,
                  bool blockkindsb = true)
        : blockStart_(0), rowEnd_(0), rowStart_(0), BlockkindSB(blockkindsb)
    {
        variableColOffsets_=new std::vector<int>();
        matrix_=new Eigen::MatrixXd();
        if (BlockkindSB)
        {
            matrix_->resize(matrix.rows(), matrix.cols());
            matrix_->triangularView<Eigen::Upper>() =
                matrix.triangularView<Eigen::Upper>();
            SfillOffsets(dimensions.begin(), dimensions.end(), appendOneDimension);
            if (matrix_->rows() != matrix_->cols())
                throw std::invalid_argument(
                    "Requested to create a SymmetricBlockMatrix from a non-square "
                    "matrix.");
            if (variableColOffsets_->back() != matrix_->cols())
                throw std::invalid_argument(
                    "Requested to create a SymmetricBlockMatrix with dimensions that "
                    "do not sum to the total size of the provided matrix.");
            assertInvariants();
        }
        else
        {
        }
    }

    SVBlockMatrix &operator=(const SVBlockMatrix &rObj)
    {
        this->BlockkindSB = rObj.BlockkindSB;
        //this->matrix_ = rObj.matrix_;
       this->matrix_->resize(rObj.matrix_->rows(), rObj.matrix_->cols());
      // this->matrix_=*rObj.matrix_;
      Eigen::MatrixXd& bb=*matrix_;
      bb=*rObj.matrix_;


        //this->variableColOffsets_ = rObj.variableColOffsets_;
        this->variableColOffsets_->clear();
        for(int keyi:*rObj.variableColOffsets_)
        {
          this->variableColOffsets_->push_back(keyi);
        }
        this->rowStart_ = rObj.rowStart_;
        this->rowEnd_ = rObj.rowEnd_;
        this->blockStart_ = rObj.blockStart_;
        return *this;
    }


    SVBlockMatrix(const SVBlockMatrix &rObj):BlockkindSB(rObj.BlockkindSB),//matrix_(rObj.matrix_),
        rowStart_(rObj.rowStart_),rowEnd_(rObj.rowEnd_),blockStart_(rObj.blockStart_)
    {
       variableColOffsets_=new std::vector<int>();
        matrix_=new Eigen::MatrixXd(rObj.matrix_->rows(),rObj.matrix_->cols());
      /*for(int i=0;i<matrix_->rows();i++)
      {
          for(int j=0;j<matrix_->cols();j++)
          {
              matrix_->operator()(i,j)=rObj.matrix_->operator()(i,j);
          }
      }*/
      Eigen::MatrixXd& bb=*matrix_;
      bb=*rObj.matrix_;



       for(int keyi:*rObj.variableColOffsets_)
        {
          this->variableColOffsets_->push_back(keyi);
        }
    }

    void setSVBlockMatrix(const SVBlockMatrix &rObj)
    {
        this->BlockkindSB = rObj.BlockkindSB;
        //this->matrix_ = rObj.matrix_;
         this->matrix_->resize(rObj.matrix_->rows(), rObj.matrix_->cols());
       //this->matrix_<<*rObj.matrix_;
       Eigen::MatrixXd& bb=*matrix_;
      bb=*rObj.matrix_;


        //this->variableColOffsets_ = rObj.variableColOffsets_;
         this->variableColOffsets_->clear();
        for(int keyi:*rObj.variableColOffsets_)
        {
          this->variableColOffsets_->push_back(keyi);
        }
        this->rowStart_ = rObj.rowStart_;
        this->rowEnd_ = rObj.rowEnd_;
        this->blockStart_ = rObj.blockStart_;
      //  return *this;
    }

    void assertInvariants() const
    {
        if (BlockkindSB)
            assert(matrix_->rows() == matrix_->cols());
        assert(matrix_->cols() == variableColOffsets_->back());
        assert(blockStart_ < (int)variableColOffsets_->size());

        assert(rowStart_ <= matrix_->rows());
        assert(rowEnd_ <= matrix_->rows());
        assert(rowStart_ <= rowEnd_);
    }

    /// Row size
    int rows() const
    {
        assertInvariants();
        if (BlockkindSB)
            return variableColOffsets_->back() - variableColOffsets_->at(blockStart_);
        else
            return rowEnd_ - rowStart_;
    }

    /// Column size
    int cols() const
    {
        if (BlockkindSB)
            return rows();
        else
            return variableColOffsets_->back() - variableColOffsets_->at(blockStart_);
    }

    /// Block count
    int nBlocks() const
    {
        if (BlockkindSB)
            return nActualBlocks() - blockStart_;
        else
            return variableColOffsets_->size() - 1 - blockStart_;
    }

    // S
    /// Number of dimensions for variable on this diagonal block.
    int SgetDim(int block) const
    {
        if (BlockkindSB)
            return calcIndices(block, block, 1, 1)[2];
        return 0;
    }

    // S
    /// @name Block getter methods.
    /// @{

    /// Get a copy of a block (anywhere in the matrix).
    /// This method makes a copy - use the methods below if performance is
    /// critical.
    Eigen::MatrixXd Sblock(int I, int J) const
    {
        if (I == J)
        {
            return SdiagonalBlock(I);
        }
        else if (I < J)
        {
            return SaboveDiagonalBlock(I, J);
        }
        else
        {
            return SaboveDiagonalBlock(J, I).transpose();
        }
    }

    // S
    /// Return the J'th diagonal block as a self adjoint view.
    Eigen::SelfAdjointView<Eigen::Block<Eigen::MatrixXd>, Eigen::Upper>
    SdiagonalBlock(int J)
    {
        return Sblock_(J, J).selfadjointView<Eigen::Upper>();
    }

    // S
    /// Return the J'th diagonal block as a self adjoint view.
    Eigen::SelfAdjointView<Eigen::Block<const Eigen::MatrixXd>, Eigen::Upper>
    SdiagonalBlock(int J) const
    {
        return Sblock_(J, J).selfadjointView<Eigen::Upper>();
    }

    // S
    /// Get the diagonal of the J'th diagonal block.
    Eigen::VectorXd diagonal(int J) const
    {
        return Sblock_(J, J).diagonal();
    }

    // S
    /// Get block above the diagonal (I, J).
    Eigen::Block<const Eigen::MatrixXd> SaboveDiagonalBlock(int I, int J) const
    {
        assert(I < J);
        return Sblock_(I, J);
    }

    // S
    /// Return the square sub-matrix that contains blocks(i:j, i:j).
    Eigen::SelfAdjointView<Eigen::Block<const Eigen::MatrixXd>, Eigen::Upper>
    selfadjointView(int I, int J) const
    {
        assert(J > I);
        return Sblock_(I, I, J - I, J - I).selfadjointView<Eigen::Upper>();
    }

    // S
    /// Return the square sub-matrix that contains blocks(i:j, i:j) as a
    /// triangular view.
    Eigen::TriangularView<Eigen::Block<const Eigen::MatrixXd>, Eigen::Upper>
    triangularView(int I, int J) const
    {
        assert(J > I);
        return Sblock_(I, I, J - I, J - I).triangularView<Eigen::Upper>();
    }

    // S
    /// Get a range [i,j) from the matrix. Indices are in block units.
    Eigen::Block<const Eigen::MatrixXd> SconstaboveDiagonalRange(
        int i_startBlock, int i_endBlock, int j_startBlock,
        int j_endBlock) const
    {
        assert(i_startBlock < j_startBlock);
        assert(i_endBlock <= j_startBlock);
        return Sblock_(i_startBlock, j_startBlock, i_endBlock - i_startBlock,
                       j_endBlock - j_startBlock);
    }

       Eigen::Block<const Eigen::MatrixXd> SaboveDiagonalRange(
        int i_startBlock, int i_endBlock, int j_startBlock,
        int j_endBlock) const
    {
        assert(i_startBlock < j_startBlock);
        assert(i_endBlock <= j_startBlock);
        return Sblock_(i_startBlock, j_startBlock, i_endBlock - i_startBlock,
                       j_endBlock - j_startBlock);
    }

    // S
    /// Get a range [i,j) from the matrix. Indices are in block units.
    Eigen::Block<Eigen::MatrixXd> SaboveDiagonalRange(int i_startBlock,
            int i_endBlock,
            int j_startBlock,
            int j_endBlock)
    {
        assert(i_startBlock < j_startBlock);
        assert(i_endBlock <= j_startBlock);
        return Sblock_(i_startBlock, j_startBlock, i_endBlock - i_startBlock,
                       j_endBlock - j_startBlock);
    }

    /// @}
    /// @name Block setter methods.
    /// @{

    // S
    /// Set a diagonal block. Only the upper triangular portion of `xpr` is
    /// evaluated.
    // template <typename XprType>
    void SsetDiagonalBlock(int I, const Eigen::MatrixXd &xpr)
    {
        Sblock_(I, I).triangularView<Eigen::Upper>() =
            xpr.triangularView<Eigen::Upper>();
    }

    // S
    /// Set an off-diagonal block. Only the upper triangular portion of `xpr` is
    /// evaluated.
    // template <typename XprType>
    void SsetOffDiagonalBlock(int I, int J, const Eigen::MatrixXd &xpr)
    {
        assert(I != J);
        if (I < J)
        {
            Sblock_(I, J) = xpr;
        }
        else
        {
            Sblock_(J, I) = xpr.transpose();
        }
    }

    // S
    /// Increment the diagonal block by the values in `xpr`. Only reads the upper
    /// triangular part of `xpr`.
    // template <typename XprType>
    void SupdateDiagonalBlock(int I, const Eigen::MatrixXd &xpr)
    {
        // TODO(gareth): Eigen won't let us add triangular or self-adjoint views
        // here, so we do it manually.
        auto dest = Sblock_(I, I);
        // cout<<"dest.rows()"<<dest.rows()<<endl<<"xpr.rows()"<<xpr.rows()<<endl;
        assert(dest.rows() == xpr.rows());
        assert(dest.cols() == xpr.cols());
        for (int col = 0; col < dest.cols(); ++col)
        {
            for (int row = 0; row <= col; ++row)
            {
                dest(row, col) += xpr(row, col);
            }
        }
    }

    // S

    /// Update an off diagonal block.
    /// NOTE(emmett): This assumes noalias().
    // template <typename XprType>
    void SupdateOffDiagonalBlock(int I, int J, const Eigen::MatrixXd &xpr)
    {
        assert(I != J);
        if (I < J)
        {
            Sblock_(I, J).noalias() += xpr;
        }
        else
        {
            Sblock_(J, I).noalias() += xpr.transpose();
        }
    }

    /// @}
    /// @name Accessing the full matrix.
    /// @{

    // S
    /// Get self adjoint view.
    Eigen::SelfAdjointView<Eigen::Block<Eigen::MatrixXd>, Eigen::Upper>
    selfadjointView()
    {
        return Sfull().selfadjointView<Eigen::Upper>();
    }

    // S
    /// Get self adjoint view.
    Eigen::SelfAdjointView<Eigen::Block<const Eigen::MatrixXd>, Eigen::Upper>
    selfadjointView() const
    {
        return Sfull().selfadjointView<Eigen::Upper>();
    }

    // S
    /// Set the entire active matrix. Only reads the upper triangular part of
    /// `xpr`.
    // template <typename XprType>
    void SsetFullMatrix(const Eigen::MatrixXd &xpr)
    {
        Sfull().triangularView<Eigen::Upper>() = xpr.triangularView<Eigen::Upper>();
    }

    // S
    /// Set the entire active matrix zero.
    void SsetZero()
    {
        Sfull().triangularView<Eigen::Upper>().setZero();
    }

    // S
    /// Negate the entire active matrix.
    void negate()
    {
        Sfull().triangularView<Eigen::Upper>() *= -1.0;
    }

    // S
    /// Invert the entire active matrix in place.
    void SinvertInPlace()
    {
        const auto identity = Eigen::MatrixXd::Identity(rows(), rows());
        Sfull().triangularView<Eigen::Upper>() =
            selfadjointView().llt().solve(identity).triangularView<Eigen::Upper>();
    }

    /// @}

    /// Retrieve or modify the first logical block, i.e. the block referenced by
    /// block index 0. Blocks before it will be inaccessible, except by accessing
    /// the underlying matrix using matrix().
    int &blockStart()
    {
        return blockStart_;
    }

    /// Retrieve the first logical block, i.e. the block referenced by block index
    /// 0. Blocks before it will be inaccessible, except by accessing the
    /// underlying matrix using matrix().
    int blockStart() const
    {
        return blockStart_;
    }

    /**
    * Given the augmented Hessian [A1'A1 A1'A2 A1'b
    *                              A2'A1 A2'A2 A2'b
    *                               b'A1  b'A2  b'b]
    * on x1 and x2, does partial Cholesky in-place to obtain [R Sd;0 L]  such
    * that R'R  = A1'A1 R'Sd = [A1'A2 A1'b] L'L is the augmented Hessian on the
    * the separator x2 R and Sd can be interpreted as a GaussianConditional |R*x1
    * + S*x2 - d]^2
    */
    void sbmcholeskyPartial(int nFrontals);

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
    // template<typename CONTAINER>
    SVBlockMatrix(std::vector<int> &dimensions, int height,
                  bool appendOneDimension, bool BlockkindSB) //variableColOffsets_(dimensions.size() + (appendOneDimension ? 2 : 1)),
        :rowStart_(0),
          rowEnd_(height),
          blockStart_(0),
          BlockkindSB(false)
    {
        // fillOffsets((dimensions.))
        //  std::vector<int>::iterator dbegin=dimensions.begin();
        //  std::vector<int>::iterator dend=dimensions.end();
        // fillOffsets(dbegin, dend, appendOneDimension);
        variableColOffsets_=new std::vector<int>(dimensions.size() + (appendOneDimension ? 2 : 1));
        matrix_=new Eigen::MatrixXd();
        VfillOffsets(dimensions.begin(), dimensions.end(), appendOneDimension);
        //  printf("variableColOffsets_.size():%d\n",variableColOffsets_.size());
        matrix_->resize(height, variableColOffsets_->back());
        assertInvariants();
    }

    // V
    /** Construct from a container of the sizes of each vertical block and a
    * pre-prepared matrix. */
    // template<typename CONTAINER, typename DERIVED>
    SVBlockMatrix(std::vector<int> &dimensions,
                  const Eigen::MatrixBase<Eigen::MatrixXd> &matrix,
                  bool appendOneDimension = false, bool BlockkindSB = false)//matrix_(matrix),//variableColOffsets_(dimensions.size() + (appendOneDimension ? 2 : 1)),
         :  rowStart_(0),
          rowEnd_(matrix.rows()),
          blockStart_(0),
          BlockkindSB(false)
    {
        variableColOffsets_=new std::vector<int>(dimensions.size() + (appendOneDimension ? 2 : 1));
        matrix_=new Eigen::MatrixXd();
        VfillOffsets(dimensions.begin(), dimensions.end(), appendOneDimension);
        if (variableColOffsets_->back() != matrix_->cols())
            throw std::invalid_argument(
                "Requested to create a VerticalBlockMatrix with dimensions that do "
                "not sum to the total columns of the provided matrix.");
        assertInvariants();
    }

    /** Construct from iterator over the sizes of each vertical block. */

    //  template<typename ITERATOR>
    SVBlockMatrix(std::vector<int>::iterator firstBlockDim,
                  std::vector<int>::iterator lastBlockDim, int height,
                  bool appendOneDimension = false, bool BlockkindSB = false)//variableColOffsets_((lastBlockDim - firstBlockDim) +(appendOneDimension ? 2 : 1)),
        :rowStart_(0),
          rowEnd_(height),
          blockStart_(0),
          BlockkindSB(false)
    {
         variableColOffsets_=new std::vector<int>((lastBlockDim - firstBlockDim) +
                              (appendOneDimension ? 2 : 1));
        VfillOffsets(firstBlockDim, lastBlockDim, appendOneDimension);
        matrix_->resize(height, variableColOffsets_->back());
        assertInvariants();
    }
    SVBlockMatrix(std::vector<int>::const_iterator firstBlockDim,
                  std::vector<int>::const_iterator lastBlockDim, int height,
                  bool appendOneDimension = false, bool BlockkindSB = false)//: variableColOffsets_((lastBlockDim - firstBlockDim) +(appendOneDimension ? 2 : 1)),
          :rowStart_(0),
          rowEnd_(height),
          blockStart_(0),
          BlockkindSB(false)
    {
    variableColOffsets_=new std::vector<int>((lastBlockDim - firstBlockDim) +
                              (appendOneDimension ? 2 : 1));
        VfillOffsets(firstBlockDim, lastBlockDim, appendOneDimension);
        matrix_->resize(height, variableColOffsets_->back());
        assertInvariants();
    }

    /** Access a single block in the underlying matrix with read/write access
    * //operator()*/
    Eigen::Block<Eigen::MatrixXd> Vblock(int blocknum)
    {
        return Vrange(blocknum, blocknum + 1);
    }

    Eigen::Block<const  Eigen::MatrixXd> Vblock(int blocknum)  const
    {
        return Vrange(blocknum, blocknum + 1);
    }

    /** Access a const block view //operator()*/
    /*Eigen::Block<const Eigen::MatrixXd> Vconstblock(int blocknum) const
    {
        return Vconstrange(blocknum, blocknum + 1);
    }*/

    /** access ranges of blocks at a time */
    Eigen::Block<Eigen::MatrixXd> Vrange(int startBlock, int endBlock)
    {
        assertInvariants();
        int actualStartBlock = startBlock + blockStart_;
        int actualEndBlock = endBlock + blockStart_;
        if (startBlock != 0 || endBlock != 0)
        {
            checkBlock(actualStartBlock);
            assert(actualEndBlock < (int)variableColOffsets_->size());
        }
        const int startCol = variableColOffsets_->at(actualStartBlock);
        const int rangeCols = variableColOffsets_->at(actualEndBlock)- startCol;
        return matrix_->block(rowStart_, startCol, this->rows(), rangeCols);
    }
    Eigen::Block<const Eigen::MatrixXd> Vrange(int startBlock,
            int endBlock) const
    {
        assertInvariants();
        int actualStartBlock = startBlock + blockStart_;
        int actualEndBlock = endBlock + blockStart_;
        if (startBlock != 0 || endBlock != 0)
        {
            checkBlock(actualStartBlock);
            assert(actualEndBlock < (int)variableColOffsets_->size());
        }
        const int startCol = variableColOffsets_->at(actualStartBlock);
        const int rangeCols = variableColOffsets_->at(actualEndBlock) - startCol;
        return ((const Eigen::MatrixXd*)(matrix_))->block(rowStart_, startCol, this->rows(), rangeCols);
    }

   /*
    Eigen::Block<const Eigen::MatrixXd> Vconstrange(int startBlock,
            int endBlock) const
    {
        assertInvariants();
        int actualStartBlock = startBlock + blockStart_;
        int actualEndBlock = endBlock + blockStart_;
        if (startBlock != 0 || endBlock != 0)
        {
            checkBlock(actualStartBlock);
            assert(actualEndBlock < (int)variableColOffsets_.size());
        }
        const int startCol = variableColOffsets_[actualStartBlock];
        const int rangeCols = variableColOffsets_[actualEndBlock] - startCol;
        return matrix_.block(rowStart_, startCol, this->rows(), rangeCols);
    }
*/
    /** Return the full matrix, *not* including any portions excluded by
    * rowStart(), rowEnd(), and firstBlock() */
    Eigen::Block<Eigen::MatrixXd> Vfull()
    {
        return Vrange(0, nBlocks());
    }

    /** Return the full matrix, *not* including any portions excluded by
    * rowStart(), rowEnd(), and firstBlock() */
    const Eigen::Block<const Eigen::MatrixXd> Vfull() const
    {
        return Vrange(0, nBlocks());
    }

    int Voffset(int block) const
    {
        assertInvariants();
        int actualBlock = block + blockStart_;
        checkBlock(actualBlock);
        return variableColOffsets_->at(actualBlock);
    }
    /** Get the apparent first row of the underlying matrix for all operations */
    const int &rowStart() const
    {
        return rowStart_;
    }

    /** Get or set the apparent first row of the underlying matrix for all
    * operations */
    int &rowStart()
    {
        return rowStart_;
    }

    /** Get the apparent last row (exclusive, i.e. rows() == rowEnd() -
    * rowStart()) of the underlying matrix for all operations */
    const int &rowEnd() const
    {
        return rowEnd_;
    }

    /** Get or set the apparent last row (exclusive, i.e. rows() == rowEnd() -
    * rowStart()) of the underlying matrix for all operations */
    int &rowEnd()
    {
        return rowEnd_;
    }

    /** Get the apparent first block for all operations */
    const int &firstBlock() const
    {
        return blockStart_;
    }

    /** Get or set the apparent first block for all operations */
    int &firstBlock()
    {
        return blockStart_;
    }

    /** Access to full matrix (*including* any portions excluded by rowStart(),
    * rowEnd(), and firstBlock()) */
    const Eigen::MatrixXd &matrix() const
    {
        return *matrix_;
    }

    /** Non-const access to full matrix (*including* any portions excluded by
    * rowStart(), rowEnd(), and firstBlock()) */
    //Eigen::MatrixXd &changematrix()
   // {
   //     return matrix_;
   // }

public:
    // S
    /// Get an offset for a block index (in the active view).
    int Soffset(int block) const
    {
        assert(block >= 0);
        const int actual_index = block + blockStart();
        assert(actual_index < nOffsets());
        return variableColOffsets_->at(actual_index);
    }

protected:
    // S
    /// Number of offsets in the full matrix.
    int nOffsets() const
    {
        return variableColOffsets_->size();
    }

    // S
    /// Number of actual blocks in the full matrix.
    int nActualBlocks() const
    {
        return nOffsets() - 1;
    }

    // S
    /// Get an arbitrary block from the matrix. Indices are in block units.
    Eigen::Block<const Eigen::MatrixXd> Sblock_(int iBlock, int jBlock,
            int blockRows = 1,
            int blockCols = 1) const
    {
        const std::array<int, 4> indices =
            calcIndices(iBlock, jBlock, blockRows, blockCols);
        return  ((const Eigen::MatrixXd*)(matrix_))->block(indices[0], indices[1], indices[2], indices[3]);
    }

    // S
    /// Get an arbitrary block from the matrix. Indices are in block units.
    Eigen::Block<Eigen::MatrixXd> Sblock_(int iBlock, int jBlock,
                                          int blockRows = 1, int blockCols = 1)
    {
        const std::array<int, 4> indices =
            calcIndices(iBlock, jBlock, blockRows, blockCols);
        return matrix_->block(indices[0], indices[1], indices[2], indices[3]);
    }

    // S
    /// Get the full matrix as a block.
    Eigen::Block<const Eigen::MatrixXd> Sfull() const
    {
        return Sblock_(0, 0, nBlocks(), nBlocks());
    }

    // S
    /// Get the full matrix as a block.
    Eigen::Block<Eigen::MatrixXd> Sfull()
    {
        return Sblock_(0, 0, nBlocks(), nBlocks());
    }

    /// Compute the indices into the underlying matrix for a given block.
    std::array<int, 4> calcIndices(int iBlock, int jBlock, int blockRows,
                                   int blockCols) const
    {
        assert(blockRows >= 0);
        assert(blockCols >= 0);

        // adjust indices to account for start and size of blocks
        const int denseI = Soffset(iBlock);
        const int denseJ = Soffset(jBlock);
        const int denseRows = Soffset(iBlock + blockRows) - denseI;
        const int denseCols = Soffset(jBlock + blockCols) - denseJ;
        return {{denseI, denseJ, denseRows, denseCols}};
    }

    // template<typename ITERATOR>
    void SfillOffsets(std::vector<int>::const_iterator firstBlockDim,
                      std::vector<int>::const_iterator lastBlockDim,
                      bool appendOneDimension)
    {
        variableColOffsets_->resize((lastBlockDim - firstBlockDim) + 1 +
                                   (appendOneDimension ? 1 : 0));
        variableColOffsets_->at(0) = 0;
        int j = 0;
        for (std::vector<int>::const_iterator dim = firstBlockDim;
                dim != lastBlockDim; ++dim)
        {
            variableColOffsets_->at(j + 1)= variableColOffsets_->at(j) + *dim;
            ++j;
        }
        if (appendOneDimension)
        {
            variableColOffsets_->at(j + 1)= variableColOffsets_->at(j) + 1;
            ++j;
        }
    }

    /// V
    void checkBlock(int block) const
    {
        static_cast<void>(block); // Disable unused varibale warnings.
        assert(matrix_->cols() == variableColOffsets_->back());
        assert(block < (int)variableColOffsets_->size() - 1);
        // cout<<matrix_.cols()<<endl;
        // cout<<variableColOffsets_[block]<<endl;
        // cout<< variableColOffsets_[block+1]<<endl;

        assert(variableColOffsets_->at(block) < matrix_->cols() &&
               variableColOffsets_->at(block + 1) <= matrix_->cols());
    }

    // template<typename ITERATOR>
    void VfillOffsets(std::vector<int>::iterator firstBlockDim,
                      std::vector<int>::iterator lastBlockDim,
                      bool appendOneDimension)
    {
        variableColOffsets_->at(0) = 0;
        int j = 0;
        for (std::vector<int>::iterator dim = firstBlockDim; dim != lastBlockDim;
                ++dim, ++j)
            variableColOffsets_->at(j + 1) = variableColOffsets_->at(j) + *dim;
        if (appendOneDimension)
            variableColOffsets_->at(j + 1) = variableColOffsets_->at(j) + 1;
    }
    void VfillOffsets(std::vector<int>::const_iterator firstBlockDim,
                      std::vector<int>::const_iterator lastBlockDim,
                      bool appendOneDimension)
    {
        variableColOffsets_->at(0) = 0;
        int j = 0;
        for (std::vector<int>::const_iterator dim = firstBlockDim; dim != lastBlockDim;
                ++dim, ++j)
            variableColOffsets_->at(j + 1) = variableColOffsets_->at(j)+ *dim;
        if (appendOneDimension)
            variableColOffsets_->at(j + 1) = variableColOffsets_->at(j) + 1;
    }
};

//
// StoS
SVBlockMatrix LikeActiveViewOfSS(const SVBlockMatrix &other);
// VtoS
SVBlockMatrix* LikeActiveViewOfVS(const SVBlockMatrix *other);
SVBlockMatrix  LikeActiveViewOfVS(const SVBlockMatrix& other);

// VtoV
/** Copy the block structure and resize the underlying matrix, but do not copy
 * the matrix data. If blockStart(), rowStart(), and/or rowEnd() have been
 * modified, this copies the structure of the corresponding matrix view. In the
 * destination VerticalBlockView, blockStart() and rowStart() will thus be 0,
 * rowEnd() will be cols() of the source VerticalBlockView, and the underlying
 * matrix will be the size of the view of the source matrix.  */
SVBlockMatrix LikeActiveViewOfVV(const SVBlockMatrix &other);

// StoV
SVBlockMatrix* LikeActiveViewOfSV(const SVBlockMatrix *rhs, int height);
SVBlockMatrix LikeActiveViewOfSV(const SVBlockMatrix& rhs, int height);
#endif // SVBLOCKMATRIX_H_INCLUDED
