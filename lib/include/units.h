#pragma once

#define EIGEN_DONT_VECTORIZE
#include <iostream>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <tbb/tbb.h>

// todo fix this stuff, make it a dropdown option
#define QUADRATIC_REGIONS
//#define AFFINE_REGIONS

#ifdef QUADRATIC_REGIONS
#define REDUCED_DOF 26
#endif
#ifdef AFFINE_REGIONS
#define REDUCED_DOF 11
#endif

using exint = int64_t;
using SolveReal = double;
using Index = exint;
//using UT_Vector3SR = UT_Vector3T<SolveReal>;
using Vector = Eigen::VectorXd;
using ReducedMatrix = Eigen::Matrix<SolveReal, REDUCED_DOF, REDUCED_DOF>;
using ColumnVector = Eigen::Matrix<SolveReal, REDUCED_DOF, 1>;
using RowVector = Eigen::Matrix<SolveReal, 1, REDUCED_DOF>;
using ReducedVector = Eigen::Matrix<SolveReal, REDUCED_DOF, 1>;
using ReducedVectorI = Eigen::Matrix<Index, REDUCED_DOF, 1>;
using SparseMatrix = Eigen::SparseMatrix<SolveReal, Eigen::RowMajor>;
using SparseMatrixColMaj = Eigen::SparseMatrix<SolveReal, Eigen::ColMajor>;
using Triplets = std::vector<Eigen::Triplet<SolveReal>>;

// compile time
static const char compileDateTime[] = { __DATE__ " / " __TIME__ };

enum AABBSide
{
	X_MIN = 0,
	Y_MIN = 1,
	Z_MIN = 2,
	X_MAX = 3,
	Y_MAX = 4,
	Z_MAX = 5
};

enum class PreconditionerType
{
    IDENTITY = 1,
    EQ_14 = 2,
    GS_SMOOTHER = 3,
	SPD_GS_SMOOTHER = 4
};

enum class MaterialLabels
{
	UNASSIGNED = -1,
	UNSOLVED = -2,
	GENERICFLUID = -3,
	ACTIVEFLUID = -4,
	SOLID = -5,
	REDUCED = -6,
	UNVISITED = -7,
	VISITED = -8,
	BOUNDARY = -9
};

namespace HDK_PolyStokes_Options
{
	enum FloatPrecision
	{
		FLOAT32,
		FLOAT64
	};

	enum MatrixScheme
	{
		PRESSURE_STRESS,
		ALL_DOFS,
		PRESSURE_VELOCITY,
		ALL_DOFS_EXPLICIT_INTERIOR_STRESS,
		PRESSURE_STRESS_EXPLICIT_INTERIOR_STRESS
	};

	enum SolverType
	{
		PCG_MATRIX_VECTOR_PRODUCTS,
		EIGEN,
		PCG_DIRECT_PRODUCTS,
		SPLIT,
		ADMM_FULL,
		DEBUG_1,
		DEBUG_2
	};
}