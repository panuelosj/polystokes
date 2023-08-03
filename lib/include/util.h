#pragma once

#include "units.h"

inline SparseMatrix
cwiseInverse(
    const SparseMatrix& matrix
)
{
    // compute component-wise inverse of a sparse matrix, ignoring all zeroes
    std::vector<Eigen::Triplet<SolveReal>> triplets;
    for (int k = 0; k < matrix.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(matrix, k); it; ++it) {
            triplets.emplace_back(it.row(), it.col(), 1./it.value());
        }
    }
    SparseMatrix retval;
    retval.resize(matrix.rows(), matrix.cols());
    retval.setFromTriplets(triplets.begin(), triplets.end());

    return retval;
}

inline SparseMatrix
extractDiagonal(
    const SparseMatrix& matrix
)
{
    // extracts only the diagonal entries of a matrix
    std::vector<Eigen::Triplet<SolveReal>> triplets;
    for (int k = 0; k < matrix.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(matrix, k); it; ++it) {
            if (it.row() == it.col())
                triplets.emplace_back(it.row(), it.col(), it.value());
        }
    }
    SparseMatrix retval;
    retval.resize(matrix.rows(), matrix.cols());
    retval.setFromTriplets(triplets.begin(), triplets.end());

    return retval;
}

inline SparseMatrix
fillEmptyDiagonalEntries(
    const SparseMatrix& matrix
)
{
    // extracts only the diagonal entries of a matrix
    std::vector<Eigen::Triplet<SolveReal>> triplets;
    for (int k = 0; k < matrix.outerSize(); ++k)
    {
        bool hasDiag = false;
        for (SparseMatrix::InnerIterator it(matrix, k); it; ++it) {
            if (it.row() == it.col()) {
                triplets.emplace_back(it.row(), it.col(), it.value());
                hasDiag = true;
                break;
            }
        }
        if (!hasDiag)
            triplets.emplace_back(k, k, 1.0);
    }
    SparseMatrix retval;
    retval.resize(matrix.rows(), matrix.cols());
    retval.setFromTriplets(triplets.begin(), triplets.end());

    return retval;
}

inline Vector
gaussSeidelIteration(
    const SparseMatrix& A, 
    const Vector b, 
    const Vector guess, 
    const exint iterations
)
{
    Vector retval = guess;
    exint n = A.rows();
    exint m = A.cols();
    assert(n == m);

    for (int q = 0; q < iterations; ++q)
    {
        for (int i = 0; i < A.outerSize(); ++i)
        {
            SolveReal sigma = 0.;
            SolveReal diag = 1.;        // todo maybe a more elegant way to avoid zero diags
            for (SparseMatrix::InnerIterator it(A, i); it; ++it) {
                if (it.row() != it.col())
                    sigma = sigma + it.value();
                else
                    diag = it.value();
            }
            retval(i) = (b(i) - sigma) / diag;
        }
    }

    return retval;
}

inline void
addBlockEntriesToTriplets(
    const SparseMatrix& matrix,
    std::vector<Eigen::Triplet<SolveReal>>& triplets,
    exint offsetRows,
    exint offsetCols,
    SolveReal scalarFactor
)
{
    for (int k = 0; k < matrix.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(matrix, k); it; ++it) {
            triplets.emplace_back(it.row() + offsetRows, it.col() + offsetCols, scalarFactor * it.value());
        }
    }
}

inline void
addBlockEntriesToTripletsTranspose(
    const SparseMatrix& matrix,
    std::vector<Eigen::Triplet<SolveReal>>& triplets,
    exint offsetRows,
    exint offsetCols,
    SolveReal scalarFactor
)
{
    for (int k = 0; k < matrix.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(matrix, k); it; ++it) {
            triplets.emplace_back(it.col() + offsetRows, it.row() + offsetCols, scalarFactor * it.value());
        }
    }
}

// is it also possible to do this by compressing the sparse matrix into a bunch of REDUCED_DOF x nActiveDofsTouched dense matrices?

//#define NEW_CODE
#define OLD_CODE

//template<typename SparseMatrixT = SparseMatrix>
inline Vector
manualMatrixTransposeVectorDistribute(
    SparseMatrix mat,
    Vector vec
)
{
    // applies mat'*vec (note the transpose)
    // where mat is row-major
    assert(mat.IsRowMajor);

    Index nDofs = mat.cols();
    Vector retval(nDofs);
    retval.setZero();
    Index nRegions = floor(vec.size() / REDUCED_DOF);

#ifdef NEW_CODE
    // todo threadcount
    Index threadCount = 2;
    Index tbbGrainSize = nRegions/(4 * threadCount);

    std::vector<Eigen::SparseVector<SolveReal>> accumulators(nRegions);
    std::vector<Eigen::SparseVector<Index>> nonzeroes(nRegions);
    
    tbb::parallel_for(tbb::blocked_range<exint>(0, nRegions, tbbGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            for (Index i = range.begin(); i != range.end(); ++i) {
                Eigen::SparseVector<SolveReal>& localAccumulator = accumulators[i];
                localAccumulator.resize(nDofs);

                for (Index j = 0; j < REDUCED_DOF; ++j) {
                    Index colnum = i * REDUCED_DOF + j;

                    localAccumulator += (Eigen::SparseVector<SolveReal>) (vec(colnum) * mat.row(colnum));
                }
            }
        });
    
    for (Eigen::SparseVector<SolveReal> accumulator : accumulators)
    {
        retval += accumulator;
    }
#endif

#ifdef OLD_CODE
    for (Index i = 0; i != nRegions; ++i) {
        for (Index j = 0; j < REDUCED_DOF; ++j) {
            Index colnum = i * REDUCED_DOF + j;
            for (SparseMatrix::InnerIterator it(mat, colnum); it; ++it) {
                retval(it.col()) = retval(it.col()) + it.value() * vec(colnum);
            }
        }
    }
#endif

    return retval;
}

// not finished
inline Vector
manualMatrixTransposeVectorDistribute2(
    SparseMatrix mat,
    Vector vec1,
    Vector vec2
)
{
    // applies mat'*vec (note the transpose)
    // where mat is row-major
    assert(mat.IsRowMajor);

    Index nDofs = mat.cols();
    Vector retval(nDofs*2);
    retval.setZero();
    Index nRegions = floor(vec1.size() / REDUCED_DOF);

    for (Index i = 0; i != nRegions; ++i) {
        for (Index j = 0; j < REDUCED_DOF; ++j) {
            Index colnum = i * REDUCED_DOF + j;
            for (SparseMatrix::InnerIterator it(mat, colnum); it; ++it) {
                retval(it.col()) = retval(it.col()) + it.value() * vec1(colnum);
                retval(it.col() + nDofs) = retval(it.col() + nDofs) + it.value() * vec2(colnum);
            }
        }
    }

    return retval;
}


template<typename SparseMatrixT = SparseMatrix>
inline void
printMatrixMemoryLayout(
    const SparseMatrixT& mat,
    const int maxEntries
)
{
    std::cout << "in memory: " << std::endl;
    int i = 0;
    for (int k = 0; k < mat.outerSize() && i < maxEntries; ++k)
    {
        for (typename SparseMatrixT::InnerIterator it(mat, k); it && i < maxEntries; ++it) {
            std::cout << it.value() << " ";
            ++i;
        }
    }
    std::cout << std::endl << std::endl;
}

template<typename SparseMatrixT = SparseMatrix>
inline Index
manualCountSparseMatrixNNZ(
    const SparseMatrixT& mat
)
{
    Index i = 0;
    for (int k = 0; k < mat.outerSize(); ++k)
    {
        for (typename SparseMatrixT::InnerIterator it(mat, k); it; ++it) {
            ++i;
        }
    }

    return i;
}

template<typename SparseMatrixT = SparseMatrix>
inline void
transposeThis(
    const SparseMatrixT& mat,
    SparseMatrixT& out
)
{
    Triplets triplets;

    for (int k = 0; k < mat.outerSize(); ++k)
    {
        for (typename SparseMatrixT::InnerIterator it(mat, k); it; ++it) {
            triplets.emplace_back(it.col(), it.row(), it.value());
        }
    }

    out.resize(mat.cols(), mat.rows());
    out.setZero();
    out.setFromTriplets(triplets.begin(), triplets.end());
}

template<typename SparseMatrixT = SparseMatrix>
inline void
transposeThisColMaj(
    const SparseMatrixT& mat,
    SparseMatrixColMaj& out
)
{
    Triplets triplets;

    for (int k = 0; k < mat.outerSize(); ++k)
    {
        for (typename SparseMatrixT::InnerIterator it(mat, k); it; ++it) {
            triplets.emplace_back(it.col(), it.row(), it.value());
        }
    }

    out.resize(mat.cols(), mat.rows());
    out.setZero();
    out.setFromTriplets(triplets.begin(), triplets.end());
}

// concatenate code from https://stackoverflow.com/questions/53617647/how-would-you-merge-sparse-matrices-to-create-a-new-sparse-matrix

template<typename ScalarT, typename StorageIndexT = int>
void sparse_stack_v(
    const Eigen::SparseMatrix<ScalarT, Eigen::ColMajor>& top,
    const Eigen::SparseMatrix<ScalarT, Eigen::ColMajor>& bottom,
    Eigen::SparseMatrix<ScalarT, Eigen::ColMajor>& stacked)
{
    assert(top.cols() == bottom.cols());
    stacked.resize(top.rows() + bottom.rows(), top.cols());
    stacked.resizeNonZeros(top.nonZeros() + bottom.nonZeros());

    StorageIndexT i = 0;

    for (StorageIndexT col = 0; col < top.cols(); col++)
    {
        stacked.outerIndexPtr()[col] = i;

        for (StorageIndexT j = top.outerIndexPtr()[col]; j < top.outerIndexPtr()[col + 1]; j++, i++)
        {
            stacked.innerIndexPtr()[i] = top.innerIndexPtr()[j];
            stacked.valuePtr()[i] = top.valuePtr()[j];
        }

        for (StorageIndexT j = bottom.outerIndexPtr()[col]; j < bottom.outerIndexPtr()[col + 1]; j++, i++)
        {
            stacked.innerIndexPtr()[i] = (StorageIndexT)top.rows() + bottom.innerIndexPtr()[j];
            stacked.valuePtr()[i] = bottom.valuePtr()[j];
        }
    }
    stacked.outerIndexPtr()[top.cols()] = i;
}


template<typename ScalarT, typename StorageIndexT = int>
void sparse_stack_h(
    const Eigen::SparseMatrix<ScalarT, Eigen::RowMajor, StorageIndexT>& left,
    const Eigen::SparseMatrix<ScalarT, Eigen::RowMajor, StorageIndexT>& right,
    Eigen::SparseMatrix<ScalarT, Eigen::RowMajor, StorageIndexT>& stacked)
{
    assert(left.row() == right.row());
    stacked.resize(left.rows(), left.cols() + right.cols());
    stacked.resizeNonZeros(left.nonZeros() + right.nonZeros());

    StorageIndexT i = 0;

    for (StorageIndexT row = 0; row < left.rows(); row++)
    {
        stacked.outerIndexPtr()[row] = i;

        for (StorageIndexT j = left.outerIndexPtr()[row]; j < left.outerIndexPtr()[row + 1]; j++, i++)
        {
            stacked.innerIndexPtr()[i] = left.innerIndexPtr()[j];
            stacked.valuePtr()[i] = left.valuePtr()[j];
        }

        for (StorageIndexT j = right.outerIndexPtr()[row]; j < right.outerIndexPtr()[row + 1]; j++, i++)
        {
            stacked.innerIndexPtr()[i] = (StorageIndexT)left.cols() + right.innerIndexPtr()[j];
            stacked.valuePtr()[i] = right.valuePtr()[j];
        }
    }
    stacked.outerIndexPtr()[left.cols()] = i;
}


template<typename ScalarT, typename StorageIndexT = int>
void sparse_stack_h(
    const Eigen::SparseMatrix<ScalarT, Eigen::ColMajor, StorageIndexT>& left,
    const Eigen::SparseMatrix<ScalarT, Eigen::ColMajor, StorageIndexT>& right,
    Eigen::SparseMatrix<ScalarT, Eigen::ColMajor>& stacked)
{
    assert(left.rows() == right.rows());

    stacked.resize(left.rows(), left.cols() + right.cols());
    stacked.resizeNonZeros(left.nonZeros() + right.nonZeros());

    std::copy(left.innerIndexPtr(), left.innerIndexPtr() + left.nonZeros(), stacked.innerIndexPtr());
    std::copy(right.innerIndexPtr(), right.innerIndexPtr() + right.nonZeros(), stacked.innerIndexPtr() + left.nonZeros());

    std::copy(left.valuePtr(), left.valuePtr() + left.nonZeros(), stacked.valuePtr());
    std::copy(right.valuePtr(), right.valuePtr() + right.nonZeros(), stacked.valuePtr() + left.nonZeros());

    std::copy(left.outerIndexPtr(), left.outerIndexPtr() + left.cols(), stacked.outerIndexPtr());//dont need the last entry of A.outerIndexPtr() -- total length is AB.cols() + 1 = A.cols() + B.cols() + 1
    std::transform(right.outerIndexPtr(), right.outerIndexPtr() + right.cols() + 1, stacked.outerIndexPtr() + left.cols(), [&](StorageIndexT i) { return i + left.nonZeros(); });
}

template<typename ScalarT, typename StorageIndexT = int>
void sparse_stack_v(
    const Eigen::SparseMatrix<ScalarT, Eigen::RowMajor, StorageIndexT>& top,
    const Eigen::SparseMatrix<ScalarT, Eigen::RowMajor, StorageIndexT>& bottom,
    Eigen::SparseMatrix<ScalarT, Eigen::RowMajor>& stacked)
{
    assert(top.cols() == bottom.cols());

    stacked.resize(top.rows() + bottom.rows(), top.cols());
    stacked.resizeNonZeros(top.nonZeros() + bottom.nonZeros());
    
    std::copy(top.innerIndexPtr(), top.innerIndexPtr() + top.nonZeros(), stacked.innerIndexPtr());
    std::copy(bottom.innerIndexPtr(), bottom.innerIndexPtr() + bottom.nonZeros(), stacked.innerIndexPtr() + top.nonZeros());

    std::copy(top.valuePtr(), top.valuePtr() + top.nonZeros(), stacked.valuePtr());
    std::copy(bottom.valuePtr(), bottom.valuePtr() + bottom.nonZeros(), stacked.valuePtr() + top.nonZeros());
    
    std::copy(top.outerIndexPtr(), top.outerIndexPtr() + top.rows(), stacked.outerIndexPtr());//dont need the last entry of A.outerIndexPtr() -- total length is AB.cols() + 1 = A.cols() + B.cols() + 1
    std::transform(bottom.outerIndexPtr(), bottom.outerIndexPtr() + bottom.rows() + 1, stacked.outerIndexPtr() + top.rows(), [&](StorageIndexT i) { return i + top.nonZeros(); });
    
}

template<typename ScalarT, typename StorageIndexT = int>
void concatenate_v(
    const Eigen::SparseMatrix<ScalarT, Eigen::RowMajor, StorageIndexT>& top,
    const Eigen::SparseMatrix<ScalarT, Eigen::RowMajor, StorageIndexT>& bottom,
    Eigen::SparseMatrix<ScalarT, Eigen::RowMajor>& stacked
)
{
    assert(top.cols() == bottom.cols());

    stacked.resize(top.rows() + bottom.rows(), top.cols());
    {
        Triplets triplets;
        triplets.reserve(top.nonZeros() + bottom.nonZeros());

        addBlockEntriesToTriplets(top, triplets, 0, 0, 1.);
        addBlockEntriesToTriplets(bottom, triplets, top.rows(), 0, 1.);

        stacked.setFromTriplets(triplets.begin(), triplets.end());
    }
}

template<typename ScalarT, typename StorageIndexT = int>
void concatenate_h(
    const Eigen::SparseMatrix<ScalarT, Eigen::RowMajor, StorageIndexT>& left,
    const Eigen::SparseMatrix<ScalarT, Eigen::RowMajor, StorageIndexT>& right,
    Eigen::SparseMatrix<ScalarT, Eigen::RowMajor>& stacked
)
{
    assert(left.rows() == right.rows());

    stacked.resize(left.rows(), left.cols() + right.cols());
    {
        Triplets triplets;
        triplets.reserve(left.nonZeros() + right.nonZeros());

        addBlockEntriesToTriplets(left, triplets, 0, 0, 1.);
        addBlockEntriesToTriplets(right, triplets, 0, left.cols(), 1.);

        stacked.setFromTriplets(triplets.begin(), triplets.end());
    }
}