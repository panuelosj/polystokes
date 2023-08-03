#include "HDK_PolyStokesSolver.h"

void
HDK_PolyStokes::Solver::constructPreconditioner()
{
    //constructPreconditionerEq14();
    //constructPreconditionerGSsmoother();
    constructPreconditionerIdentity();
}

void
HDK_PolyStokes::Solver::constructPreconditionerGSsmoother()
{
    assembleReducedMassBlock();         // need Mr
    assembleReducedInvertedBlock();     // need Binv
    assembleVMatrices();

    presolver = new Preconditioner(PreconditionerType::GS_SMOOTHER);
    presolver->setupGSsmoother(
        Mc_Matrix,
        Mr_Matrix,
        Inv_Mr_plus_2JDtuDJ_Matrix,
        V_Matrix,
        G_Matrix,
        VJt_Matrix,
        JG_Matrix,
        dt
    );
}

void
HDK_PolyStokes::Solver::constructPreconditionerIdentity()
{
    presolver = new Preconditioner(PreconditionerType::IDENTITY);
}

void
HDK_PolyStokes::Solver::constructPreconditionerDiagonal()
{
    // todo
}

void
HDK_PolyStokes::Solver::constructPreconditionerEq14()
{
    // we need the reduced mass block
    assembleReducedMassBlock();
    assembleReducedMassBlockInverse();

    SparseMatrix A1;
    A1.resize(nPressures, nActiveVs + nReducedVs);
    {
        Triplets triplets;
        exint sparseNonzeroElems =
            G_Matrix.nonZeros()
            + JG_Matrix.nonZeros();
        triplets.reserve(sparseNonzeroElems);

        addBlockEntriesToTripletsTranspose(G_Matrix, triplets, 0, 0, 1.);
        addBlockEntriesToTripletsTranspose(JG_Matrix, triplets, 0, nActiveVs, 1.);

        A1.setFromTriplets(triplets.begin(), triplets.end());
    }
    SparseMatrix Dtilde;
    Dtilde.resize(nActiveVs + nReducedVs, nActiveVs + nReducedVs);
    {
        Triplets triplets;
        exint sparseNonzeroElems =
            Mc_Matrix.nonZeros()
            + Mr_Matrix.nonZeros();
        triplets.reserve(sparseNonzeroElems);

        addBlockEntriesToTriplets(Mc_Matrix, triplets, 0, 0, invDt);
        addBlockEntriesToTriplets(Mr_Matrix, triplets, nActiveVs, nActiveVs, invDt);

        Dtilde.setFromTriplets(triplets.begin(), triplets.end());
    }
    SparseMatrix DtildeInv;
    DtildeInv.resize(nActiveVs + nReducedVs, nActiveVs + nReducedVs);
    {
        Triplets triplets;
        exint sparseNonzeroElems =
            McInv_Matrix.nonZeros()
            + MrInv_Matrix.nonZeros();
        triplets.reserve(sparseNonzeroElems);

        addBlockEntriesToTriplets(McInv_Matrix, triplets, 0, 0, dt);
        addBlockEntriesToTriplets(MrInv_Matrix, triplets, nActiveVs, nActiveVs, dt);

        DtildeInv.setFromTriplets(triplets.begin(), triplets.end());
    }

    presolver = new Preconditioner(PreconditionerType::EQ_14);
    presolver->setupEq14Inv(A1, Dtilde, DtildeInv);
}

void
HDK_PolyStokes::Solver::constructPreconditionerEq6()
{
    // todo
}
