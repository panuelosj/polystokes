#include "Preconditioner.h"
// saving out matrix
#include <unsupported/Eigen/SparseExtra>

// constructor and destructor
Preconditioner::Preconditioner(
	PreconditionerType _type
) :
	myType(_type)
{
}

Preconditioner::~Preconditioner()
{
}

// performs the action of the inverse of the preconditioner
Vector Preconditioner::solve(Vector b)
{
	if (myType == PreconditionerType::EQ_14)
		return solveEq14(b);
	else if (myType == PreconditionerType::GS_SMOOTHER)
		return solveGSsmoother(b);
	else if (myType == PreconditionerType::SPD_GS_SMOOTHER)
		return solveSPDgsSmoother(b);
	else
		return solveIdentity(b);
}

void Preconditioner::setupGSsmoother(
	SparseMatrix _Mc,
	SparseMatrix _Mr,
	SparseMatrix _Binv,
	SparseMatrix _V,
	SparseMatrix _G,
	SparseMatrix _VJt,
	SparseMatrix _JG,
	SolveReal _dt
)
{
	Mc = _Mc;
	Mr = _Mr;
	Binv = _Binv;
	V = _V;
	G = _G;
	VJt = _VJt;
	JG = _JG;
	dt = _dt;

	GSsmootherReady = true;
}

Vector Preconditioner::solveGSsmoother(const Vector b)
{
	if (!GSsmootherReady) return b;
	
	Index n_us = Mc.rows();
	Index n_vs = Mr.rows();
	Index n_ps = G.cols();

	Vector z_u(n_us), z_v(n_vs), z_p(n_ps);
	Vector r_u(n_us), r_v(n_vs), r_p(n_ps);
	Vector z_up(n_us + n_ps);
	r_u = b.head(n_us);
	r_v = b.segment(n_us, n_vs);
	r_p = b.tail(n_ps);

	// step 1
	z_u = stepGSsmootherUniform(r_u, z_v, z_p);
	//z_up = stepGSsmootherUniformExact(r_u, z_v, z_p);
	//z_u = z_up.head(n_us);
	//z_p = z_up.tail(n_ps);
	// step 2
	z_v = stepGSsmootherReduced(r_v, z_u, z_p);
	// step 3
	z_u = stepGSsmootherUniform(r_u, z_v, z_p);
	//z_up = stepGSsmootherUniformExact(r_u, z_v, z_p);
	//z_u = z_up.head(n_us);
	//z_p = z_up.tail(n_ps);

	Vector retval = Vector::Zero(b.size());
	retval << z_u, z_v, z_p;

	return retval;
}

Vector Preconditioner::stepGSsmootherUniform(Vector r_u, Vector z_v, Vector z_p)
{
	// freeze v_R
	// do we also need to freeze p's since they're zero diagonal??
	// maybe we can just use 

	Vector rhs = (1./dt)*Mc*r_u + VJt*z_v - G*z_p;
	SparseMatrix mat = (1./dt)*Mc - V;

	Vector z_u = gaussSeidelIteration(mat, rhs, r_u, 16);

	return z_u;
}

Vector Preconditioner::stepGSsmootherUniformExact(Vector r_u, Vector z_v, Vector z_p)
{
	//std::cout << "here1" << std::endl;
	Index n_us = r_u.size();
	Index n_ps = z_p.size();
	Index usOffset = 0;
	Index psOffset = n_us;

	//std::cout << "here2" << std::endl;
	// first setup matrix
	SparseMatrix mat;
	mat.resize(n_us+n_ps, n_us+n_ps);
	{
		Triplets triplets;

		addBlockEntriesToTriplets(Mc, triplets, usOffset, usOffset, 1./dt);
		addBlockEntriesToTriplets(V, triplets, usOffset, usOffset, -1.);
		addBlockEntriesToTriplets(G, triplets, usOffset, psOffset, 1.);
		addBlockEntriesToTripletsTranspose(G, triplets, psOffset, usOffset, 1.);

		mat.setFromTriplets(triplets.begin(), triplets.end());
	}
	//std::cout << "here3" << std::endl;
	// setup rhs
	Vector rhs(n_us+n_ps);
	{
		Vector rhs_us(n_us);
		Vector rhs_ps(n_ps);

		rhs_us = (1./dt) * Mc * r_u + VJt * z_v;
		rhs_ps = - JG.transpose() * z_v;

		for (Index i=0; i<n_us; ++i)
			rhs(usOffset + i) = rhs_us(i);
		for (Index i=0; i<n_ps; ++i)
			rhs(psOffset + i) = rhs_ps(i);
	}
	//std::cout << "here4" << std::endl;
	// solution
	Vector solutionVector;

	// now solve
	Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper> solver;
	solver.compute(mat);
	solver.setMaxIterations(25);
	solver.setTolerance(1e-2);
	solutionVector = solver.solve(rhs);
	//std::cout << "here5" << std::endl;

	return solutionVector;
}

Vector Preconditioner::stepGSsmootherReduced(Vector r_v, Vector z_u, Vector z_p)
{
	Vector z_v = Binv * ((1./dt)* Mr * r_v + VJt.transpose()*z_u - JG*z_p);

	return z_v;
}

void Preconditioner::setupSPDgsSmoother(
	SparseMatrix _Mc,
	SparseMatrix _B,
	SparseMatrix _G,
	SparseMatrix _JG,
	SparseMatrix _Dt,
	SparseMatrix _JDt
)
{
	Mc = _Mc;
	B = _B;
	Dt = _Dt;
	JDt = _JDt;
	G = _G;
	JG = _JG;

	SparseMatrix Gt = G.transpose();
	SparseMatrix D = Dt.transpose();
	SparseMatrix GtJt = JG.transpose();
	SparseMatrix DJt = JDt.transpose();

	concatenate_h(G, Dt, cat_G_Dt);
	concatenate_v(Gt, D, cat_Gt_D);
	concatenate_h(JG, JDt, cat_JG_JDt);
	concatenate_v(GtJt, DJt, cat_GtJt_DJt);

	SPDgsSmootherReady = true;

	return;
}

Vector Preconditioner::solveSPDgsSmoother(Vector b)
{
	if (!SPDgsSmootherReady) return b;

	//Vector step1 = cat_Gt_D * Mc * cat_G_Dt * b;
	Vector step2 = -(1. / dt) * cat_GtJt_DJt * B * cat_JG_JDt * b;
	//Vector step3 = -(1./dt) * cat_Gt_D * Mc * cat_G_Dt * step2;

	return step2;
}

void Preconditioner::setupEq14Inv(
	SparseMatrix A, 
	SparseMatrix Dtilde, 
	SparseMatrix DtildeInv )
{
	Index n = A.cols();
	Index m = A.rows();

	SparseMatrix In, Im;
	In.resize(n, n);
	In.setIdentity();
	Im.resize(m, m);
	Im.setIdentity();

	SparseMatrix ADinv = A * DtildeInv;
	SparseMatrix ADinvAt = A * DtildeInv * A.transpose();
	//SparseMatrix DiagbADinvAt = SparseMatrix(ADinvAt.diagonal().asDiagonal());
	SparseMatrix DiagbADinvAt = fillEmptyDiagonalEntries(extractDiagonal(ADinvAt));
	SparseMatrix DiagbADinvAtInv = cwiseInverse(DiagbADinvAt);

	// todo check if this inverse is correct
	Eigen::saveMarket(DiagbADinvAt, "output_data/Mat_Diag.mtx");
	Eigen::saveMarket(DiagbADinvAtInv, "output_data/Mat_DiagInv.mtx");

	M1inv.resize(n + m, n + m);
	{
		std::vector<Eigen::Triplet<SolveReal>> triplets;
		exint sparseNonzeroElems =
			In.nonZeros()
			+ Im.nonZeros()
			+ ADinv.nonZeros();
		triplets.reserve(sparseNonzeroElems);

		addBlockEntriesToTriplets(In, triplets, 0, 0, 1.);
		addBlockEntriesToTriplets(Im, triplets, n, n, 1.);
		addBlockEntriesToTriplets(ADinv, triplets, n, 0, -1.);

		M1inv.setFromTriplets(triplets.begin(), triplets.end());
	}
	M2inv.resize(n + m, n + m);
	{
		std::vector<Eigen::Triplet<SolveReal>> triplets;
		exint sparseNonzeroElems =
			DtildeInv.nonZeros()
			+ DiagbADinvAtInv.nonZeros();
		triplets.reserve(sparseNonzeroElems);

		addBlockEntriesToTriplets(DtildeInv, triplets, 0, 0, 1.);
		addBlockEntriesToTriplets(DiagbADinvAtInv, triplets, n, n, -1.);

		M2inv.setFromTriplets(triplets.begin(), triplets.end());
	}
	M3inv = M1inv.transpose();

	Eigen::saveMarket(M1inv, "output_data/Mat_M1inv.mtx");
	Eigen::saveMarket(M2inv, "output_data/Mat_M2inv.mtx");
	Eigen::saveMarket(M3inv, "output_data/Mat_M3inv.mtx");

	Eq14Ready = true;
}

Vector Preconditioner::solveEq14(Vector b)
{
	if (Eq14Ready)
		return M3inv * M2inv * M1inv * b;
	else
		return b;
}

Vector Preconditioner::solveIdentity(Vector b)
{
	return b;
}