#pragma once

#include "units.h"
#include "util.h"

// saving out matrix
#include <unsupported/Eigen/SparseExtra>

template<typename SolveRealT = SolveReal>
class ApplyPressureStressMatrix
{
public:
	using SolveReal = SolveRealT;
	using SparseMatrix = Eigen::SparseMatrix<SolveRealT, Eigen::RowMajor>;
	using Vector = Eigen::Matrix<SolveRealT, -1, 1>;

	// Constructor and Destructor
	explicit ApplyPressureStressMatrix() {}
	virtual ~ApplyPressureStressMatrix() {}

	void setupDirect(SparseMatrix _A) {
		A = _A;
	}
	void setupMatrixVectorProducts(
		SolveReal _dt,
		SolveReal _invDt,
		SparseMatrix _McInv_Matrix,
		SparseMatrix _BInv_Matrix,
		SparseMatrix _uInv_Matrix,
		SparseMatrix _G_Matrix,
		SparseMatrix _JG_Matrix,
		SparseMatrix _Dt_Matrix,
		SparseMatrix _JDt_Matrix
	)
	{
		dt = _dt;
		invDt = _invDt;
		McInv_Matrix = Eigen::DiagonalMatrix<SolveReal, Eigen::Dynamic>(Eigen::Diagonal(_McInv_Matrix));
		BInv_Matrix = _BInv_Matrix;
		uInv_Matrix = Eigen::DiagonalMatrix<SolveReal, Eigen::Dynamic>(Eigen::Diagonal(_uInv_Matrix));
		G_Matrix = _G_Matrix;
		JG_Matrix = _JG_Matrix;
		Dt_Matrix = _Dt_Matrix;
		JDt_Matrix = _JDt_Matrix;

		Gt_Matrix = G_Matrix.transpose().eval();
		GtJt_Matrix = JG_Matrix.transpose().eval();
		D_Matrix = Dt_Matrix.transpose().eval();
		DJt_Matrix = JDt_Matrix.transpose().eval();

		BInv_Matrix.makeCompressed();
		G_Matrix.makeCompressed();
		JG_Matrix.makeCompressed();
		Dt_Matrix.makeCompressed();
		JDt_Matrix.makeCompressed();
		Gt_Matrix.makeCompressed();
		GtJt_Matrix.makeCompressed();
		D_Matrix.makeCompressed();
		DJt_Matrix.makeCompressed();

		concatenate_h(G_Matrix, Dt_Matrix, cat_G_Dt);
		concatenate_v(Gt_Matrix, D_Matrix, cat_Gt_D);
		concatenate_h(JG_Matrix, JDt_Matrix, cat_JG_JDt);
		concatenate_v(GtJt_Matrix, DJt_Matrix, cat_GtJt_DJt);

		n_ps = G_Matrix.cols();
		n_ts = Dt_Matrix.cols();
	}

	Vector applyDirect(Vector x)
	{
		return A * x;
	}

	Vector applyMatrixVectorProductsAlt(Vector x)
	{

		Vector ndtx = -dt * x;

		Vector V3_1 = - cat_JG_JDt * x;
		Vector V3_2 = BInv_Matrix * V3_1;
		//Vector V3 = manualMatrixTransposeVectorDistribute(cat_GtJt_DJt, V3_2);
		Vector V3 = cat_GtJt_DJt * V3_2;
		//Vector V3(V3_2.size());
		//V3 = manualMatrixTransposeVectorDistribute(cat_GtJt_DJt, V3_2);
		
		Vector V2_1 = cat_G_Dt * ndtx;
		Vector V2_2 = McInv_Matrix * V2_1;
		Vector V2 = cat_Gt_D * V2_2;

		Vector V1_1 = -0.5 * uInv_Matrix * x.tail(n_ts);
		Vector V1(V2.size());
		V1.setZero();
		V1.tail(n_ts) = V1_1;

		Vector retval(V1.size());
		retval = V1 + V2 + V3;

		return retval;
	}
	
	Vector applyMatrixVectorProducts(Vector x)
	{
		Vector x_ps = x.head(n_ps);
		Vector x_ts = x.tail(n_ts);

		// why is this so slow. this should have
		// 3 applications of G, D, JG, JDt each
		// 2 applications of Binv, McInv
		// 1 application of uInv
		// 	   233532+15000+10957 = 259 489
		// but A has 1 685 347 nnz
		// nnzs: G = 9191, D = 23751, JG = 11024, JDt = 33878
		// Binv = 1352, McInv = 6148
		// uInv = 10957
		// 
		// okay these flops actually make sense,
		// the sparse matrix - vector multiply requires nnz + 3 * number of rows for the accumulator on each row (2 +'s and 1 *)
		// it just looks really bad here since G is fairly sparse but has a large number of rows (each row only has an average of 1.5 nnz)
		Vector A11_1, A11_2, A21_1, A21_2, A12_1, A12_2, A22_1, A22_2;

#pragma omp parallel sections
		{
#pragma omp section
			{
		SparseMatrix McInv_G = McInv_Matrix * G_Matrix;
		Vector McInv_G_xps = McInv_G * x_ps;
		A11_1 = McInv_G_xps;
		A11_1 = Gt_Matrix * A11_1;	// ~2.7x
		A11_1 = -dt * A11_1;
		A21_1 = McInv_G_xps;
		A21_1 = D_Matrix * A21_1;
		A21_1 = -dt * A21_1;
			}

#pragma omp section
			{
		Vector BInv_JDt_xts = JDt_Matrix * x_ts;
		BInv_JDt_xts = BInv_Matrix * BInv_JDt_xts;
		Vector BInv_JG_xps = JG_Matrix * x_ps;
		BInv_JG_xps = BInv_Matrix * BInv_JG_xps;

		Vector tmp = -manualMatrixTransposeVectorDistribute2(JG_Matrix, BInv_JG_xps, BInv_JDt_xts);
		Index nDofs = tmp.size() / 2;
		A11_2 = tmp.head(nDofs);
		A12_2 = tmp.tail(nDofs);

		tmp = -manualMatrixTransposeVectorDistribute2(JDt_Matrix, BInv_JG_xps, BInv_JDt_xts);
		nDofs = tmp.size() / 2;
		A21_2 = tmp.head(nDofs);
		A22_2 = tmp.tail(nDofs);
			}

#pragma omp section
			{
		SparseMatrix McInv_Dt = McInv_Matrix * Dt_Matrix;
		Vector McInv_Dt_xts = McInv_Dt * x_ts;
		A12_1 = McInv_Dt_xts;
		A12_1 = -dt * Gt_Matrix * A12_1;
		A22_1 = McInv_Dt_xts;
		A22_1 = -dt * D_Matrix * A22_1;
			}

		}

		Vector A11 = A11_1 + A11_2;
		Vector A21 = A21_1 + A21_2;

		Vector A12 = A12_1 + A12_2;
		Vector A22_3 = -0.5 * uInv_Matrix * x_ts;
		Vector A22 = A22_1 + A22_2 + A22_3;

		Vector ret_ps = A11 + A12;
		Vector ret_ts = A21 + A22;
		Vector retval(n_ps + n_ts);
		retval << ret_ps, ret_ts;

		return retval;
	}
	

	Vector apply(Vector x) {
		return applyMatrixVectorProducts(x);
	}

private:
	Index n_ps, n_ts;

	// for direct
	SparseMatrix A;

	// for matrix-vector products
	Eigen::DiagonalMatrix<SolveReal, Eigen::Dynamic> McInv_Matrix, uInv_Matrix;
	SparseMatrix BInv_Matrix;
	SparseMatrix G_Matrix, JG_Matrix, Dt_Matrix, JDt_Matrix;
	SparseMatrix Gt_Matrix, GtJt_Matrix, D_Matrix, DJt_Matrix;
	SolveReal dt, invDt;

	// for alt matrix-vector products
	SparseMatrix cat_G_Dt, cat_Gt_D, cat_JG_JDt, cat_GtJt_DJt;
};