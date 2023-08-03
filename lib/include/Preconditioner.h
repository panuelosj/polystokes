#pragma once

#include "units.h"
#include "util.h"

class Preconditioner
{
public:
	// Constructor and Destructor
	explicit Preconditioner(
		PreconditionerType _type
	);
	virtual ~Preconditioner();

	Vector solve(Vector b);

	Vector solveIdentity(Vector b);

	void setupEq14Inv(
		SparseMatrix A, 
		SparseMatrix Dtilde, 
		SparseMatrix DtildeInv);
	Vector solveEq14(Vector b);

	void setupGSsmoother(
		SparseMatrix _Mc,
		SparseMatrix _Mr,
		SparseMatrix _Binv,
		SparseMatrix _V,
		SparseMatrix _G,
		SparseMatrix _VJt,
		SparseMatrix _JG,
		SolveReal _dt
	);
	Vector solveGSsmoother(Vector b);
	Vector stepGSsmootherUniform(Vector r_u, Vector z_v, Vector z_p);
	Vector stepGSsmootherUniformExact(Vector r_u, Vector z_v, Vector z_p);
	Vector stepGSsmootherReduced(Vector r_v, Vector z_u, Vector z_p);

	void setupSPDgsSmoother(
		SparseMatrix _Mc,
		SparseMatrix _B,
		SparseMatrix _G,
		SparseMatrix _JG,
		SparseMatrix _Dt, 
		SparseMatrix _JDt
	);
	Vector solveSPDgsSmoother(Vector b);

private:
	PreconditionerType myType = PreconditionerType::IDENTITY;

	// for GSsmoother
	bool GSsmootherReady = false;
	SolveReal dt;
	SparseMatrix Mc, Mr, Binv, V, G, VJt, JG;
	Vector r_vr, r_u;
	
	// for Eq 14
	bool Eq14Ready = false;
	SparseMatrix M1inv, M2inv, M3inv;

	// for SPDgsSmoother
	bool SPDgsSmootherReady = false;
	SparseMatrix B, Dt, JDt, cat_G_Dt, cat_Gt_D, cat_JG_JDt, cat_GtJt_DJt;
};

