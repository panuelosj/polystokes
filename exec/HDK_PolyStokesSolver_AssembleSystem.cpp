#include "HDK_PolyStokesSolver.h"

void
HDK_PolyStokes::Solver::assemble()
{
    // we only need to assemble if we're using a general solver
    if (mySolverType == HDK_PolyStokes_Options::SolverType::EIGEN || mySolverType == HDK_PolyStokes_Options::SolverType::PCG_DIRECT_PRODUCTS)
    {
        switch (myMatrixScheme)
        {
        case HDK_PolyStokes_Options::MatrixScheme::ALL_DOFS:
            assembleSystem();
            break;
        case HDK_PolyStokes_Options::MatrixScheme::PRESSURE_VELOCITY:
            assembleSystemVelocityPressure();
            break;
        case HDK_PolyStokes_Options::MatrixScheme::PRESSURE_STRESS:
            assembleSystemPressureStress();
            break;
        case HDK_PolyStokes_Options::MatrixScheme::ALL_DOFS_EXPLICIT_INTERIOR_STRESS:
            assembleSystemExplicitInternalStresses();
            break;
        }
    }
    else
    {
        //std::cout << "custom solver does not need matrix to be assembled" << std::endl;
        //std::cout << "we still should assemble blocks now" << std::endl;
        
        switch (mySolverType)
        {
        case HDK_PolyStokes_Options::SolverType::PCG_MATRIX_VECTOR_PRODUCTS:
            assembleSystemPressureStressFactored();
            break;
        }
    }
}

void
HDK_PolyStokes::Solver::assembleSystem()
{
    // create offsets
    exint activeVsOffset = 0;
    exint reducedVsOffset = nActiveVs;
    exint pressureOffset = nActiveVs + nReducedVs;
    exint stressOffset = nActiveVs + nReducedVs + nPressures;
    nSystemSize = nActiveVs + nReducedVs + nPressures + nStresses;

    // first assemble the reduced matrices and vectors
    //assembleReducedMassBlock();
    //assembleReducedViscosityBlock();
    assembleReducedCombinedBlock();
    assembleReducedRHSVector();

    A.resize(nSystemSize, nSystemSize);
    // now assemble the entire thing together
    {
        Triplets triplets;
        exint sparseNonzeroElems =
            Mc_Matrix.nonZeros()
            //+ Mr_plus_2JDtuDJ_Matrix.nonZeros()
            //+ Mr_Matrix.nonZeros()
            + myInteriorRegionCount * REDUCED_DOF * REDUCED_DOF
            + 2 * G_Matrix.nonZeros()
            + 2 * Dt_Matrix.nonZeros()
            + 2 * JG_Matrix.nonZeros()
            + 2 * JDt_Matrix.nonZeros()
            + uInv_Matrix.nonZeros();
        triplets.reserve(sparseNonzeroElems);

        addBlockEntriesToTriplets(Mc_Matrix, triplets, activeVsOffset, activeVsOffset, invDt);
        addBlockEntriesToTriplets(G_Matrix, triplets, activeVsOffset, pressureOffset, 1.0);
        addBlockEntriesToTriplets(Dt_Matrix, triplets, activeVsOffset, stressOffset, 1.0);

        addBlockEntriesToTriplets(Mr_plus_2JDtuDJ_Matrix, triplets, reducedVsOffset, reducedVsOffset, 1.0);
        addBlockEntriesToTriplets(JG_Matrix, triplets, reducedVsOffset, pressureOffset, 1.0);
        addBlockEntriesToTriplets(JDt_Matrix, triplets, reducedVsOffset, stressOffset, 1.0);

        addBlockEntriesToTripletsTranspose(G_Matrix, triplets, pressureOffset, activeVsOffset, 1.0);
        addBlockEntriesToTripletsTranspose(JG_Matrix, triplets, pressureOffset, reducedVsOffset, 1.0);

        addBlockEntriesToTripletsTranspose(Dt_Matrix, triplets, stressOffset, activeVsOffset, 1.0);
        addBlockEntriesToTripletsTranspose(JDt_Matrix, triplets, stressOffset, reducedVsOffset, 1.0);
        addBlockEntriesToTriplets(uInv_Matrix, triplets, stressOffset, stressOffset, -0.5);

        A.setFromTriplets(triplets.begin(), triplets.end());
        //A.makeCompressed();
    }

    b = Vector::Zero(nSystemSize);
    {
        for (exint i = 0; i < nActiveVs; ++i)
            b(activeVsOffset + i) = invDt * activeRHSVector(i);
        for (exint i = 0; i < nReducedVs; ++i)
            b(reducedVsOffset + i) = invDt * reducedRHSVector(i);
        for (exint i = 0; i < nPressures; ++i)
            b(pressureOffset + i) = pressureRHSVector(i);
        for (exint i = 0; i < nStresses; ++i)
            b(stressOffset + i) = stressRHSVector(i);
    }

    guessVector = Vector::Zero(nSystemSize);
    {
        for (exint i = 0; i < nActiveVs; ++i)
            guessVector(activeVsOffset + i) = activeGuessVector(i);
        for (exint i = 0; i < nPressures; ++i)
            guessVector(pressureOffset + i) = pressureGuessVector(i);
        for (exint i = 0; i < nStresses; ++i)
            guessVector(stressOffset + i) = stressGuessVector(i);
        for (exint i = 0; i < nReducedVs; ++i)
            guessVector(reducedVsOffset + i) = reducedGuessVector(i);
    }

    solutionVector = Vector::Zero(nSystemSize);
}

void
HDK_PolyStokes::Solver::assembleSystemAlt()
{
    // create offsets
    exint activeVsOffset = 0;
    exint reducedVsOffset = nActiveVs;
    exint pressureOffset = nActiveVs + nReducedVs;
    exint stressOffset = nActiveVs + nReducedVs + nPressures;
    nSystemSize = nActiveVs + nReducedVs + nPressures + nStresses;

    // first assemble the reduced matrices and vectors
    assembleReducedMassBlock();
    //assembleReducedViscosityBlock();
    //assembleReducedCombinedBlock();
    assembleReducedRHSVector();

    A.resize(nSystemSize, nSystemSize);
    // now assemble the entire thing together
    {
        Triplets triplets;
        exint sparseNonzeroElems =
            Mc_Matrix.nonZeros()
            //+ Mr_plus_2JDtuDJ_Matrix.nonZeros()
            //+ Mr_Matrix.nonZeros()
            + myInteriorRegionCount * REDUCED_DOF * REDUCED_DOF
            + 2 * G_Matrix.nonZeros()
            + 2 * Dt_Matrix.nonZeros()
            + 2 * JG_Matrix.nonZeros()
            + 2 * JDt_Matrix.nonZeros()
            + uInv_Matrix.nonZeros();
        triplets.reserve(sparseNonzeroElems);

        SparseMatrix JDtuDJt_Red_Matrix = 2. * JDtRed_Matrix * uRed_Matrix * JDtRed_Matrix.transpose();

        addBlockEntriesToTriplets(Mc_Matrix, triplets, activeVsOffset, activeVsOffset, invDt);
        addBlockEntriesToTriplets(G_Matrix, triplets, activeVsOffset, pressureOffset, 1.0);
        addBlockEntriesToTriplets(Dt_Matrix, triplets, activeVsOffset, stressOffset, 1.0);

        addBlockEntriesToTriplets(Mr_Matrix, triplets, reducedVsOffset, reducedVsOffset, invDt);
        //addBlockEntriesToTriplets(JDtuDJt_Red_Matrix, triplets, reducedVsOffset, reducedVsOffset, 1.0);
        addBlockEntriesToTriplets(JG_Matrix, triplets, reducedVsOffset, pressureOffset, 1.0);
        addBlockEntriesToTriplets(JDt_Matrix, triplets, reducedVsOffset, stressOffset, 1.0);

        addBlockEntriesToTripletsTranspose(G_Matrix, triplets, pressureOffset, activeVsOffset, 1.0);
        addBlockEntriesToTripletsTranspose(JG_Matrix, triplets, pressureOffset, reducedVsOffset, 1.0);

        addBlockEntriesToTripletsTranspose(Dt_Matrix, triplets, stressOffset, activeVsOffset, 1.0);
        addBlockEntriesToTripletsTranspose(JDt_Matrix, triplets, stressOffset, reducedVsOffset, 1.0);
        addBlockEntriesToTriplets(uInv_Matrix, triplets, stressOffset, stressOffset, -0.5);

        A.setFromTriplets(triplets.begin(), triplets.end());
        //A.makeCompressed();
    }

    b = Vector::Zero(nSystemSize);
    {
        for (exint i = 0; i < nActiveVs; ++i)
            b(activeVsOffset + i) = invDt * activeRHSVector(i);
        for (exint i = 0; i < nReducedVs; ++i)
            b(reducedVsOffset + i) = invDt * reducedRHSVector(i);
        for (exint i = 0; i < nPressures; ++i)
            b(pressureOffset + i) = pressureRHSVector(i);
        for (exint i = 0; i < nStresses; ++i)
            b(stressOffset + i) = stressRHSVector(i);
    }

    guessVector = Vector::Zero(nSystemSize);
    {
        for (exint i = 0; i < nActiveVs; ++i)
            guessVector(activeVsOffset + i) = activeGuessVector(i);
        for (exint i = 0; i < nPressures; ++i)
            guessVector(pressureOffset + i) = pressureGuessVector(i);
        for (exint i = 0; i < nStresses; ++i)
            guessVector(stressOffset + i) = stressGuessVector(i);
        for (exint i = 0; i < nReducedVs; ++i)
            guessVector(reducedVsOffset + i) = reducedGuessVector(i);
    }

    solutionVector = Vector::Zero(nSystemSize);
}

void
HDK_PolyStokes::Solver::assembleSystemExplicitInternalStresses()
{
    // create offsets
    exint activeVsOffset = 0;
    exint reducedVsOffset = nActiveVs;
    exint pressureOffset = nActiveVs + nReducedVs;
    exint stressOffset = nActiveVs + nReducedVs + nPressures;
    exint reducedStressOffset = nActiveVs + nReducedVs + nPressures + nStresses;
    nSystemSize = nActiveVs + nReducedVs + nPressures + nStresses + nReducedStresses;

    // first assemble the reduced matrices and vectors
    assembleReducedMassBlock();
    //assembleReducedViscosityBlock();
    //assembleReducedCombinedBlock();
    assembleReducedRHSVector();

    A.resize(nSystemSize, nSystemSize);
    // now assemble the entire thing together
    {
        Triplets triplets;
        exint sparseNonzeroElems =
            Mc_Matrix.nonZeros()
            + myInteriorRegionCount * REDUCED_DOF * REDUCED_DOF
            + 2 * G_Matrix.nonZeros()
            + 2 * Dt_Matrix.nonZeros()
            + 2 * JG_Matrix.nonZeros()
            + 2 * JDt_Matrix.nonZeros()
            + 2 * JDtRed_Matrix.nonZeros()
            + uInv_Matrix.nonZeros()
            + uInvRed_Matrix.nonZeros();
        triplets.reserve(sparseNonzeroElems);

        addBlockEntriesToTriplets(Mc_Matrix, triplets, activeVsOffset, activeVsOffset, invDt);
        addBlockEntriesToTriplets(G_Matrix, triplets, activeVsOffset, pressureOffset, 1.0);
        addBlockEntriesToTriplets(Dt_Matrix, triplets, activeVsOffset, stressOffset, 1.0);

        addBlockEntriesToTriplets(Mr_Matrix, triplets, reducedVsOffset, reducedVsOffset, 1.0);
        addBlockEntriesToTriplets(JG_Matrix, triplets, reducedVsOffset, pressureOffset, 1.0);
        addBlockEntriesToTriplets(JDt_Matrix, triplets, reducedVsOffset, stressOffset, 1.0);
        addBlockEntriesToTriplets(JDtRed_Matrix, triplets, reducedVsOffset, reducedStressOffset, 1.0);

        addBlockEntriesToTripletsTranspose(G_Matrix, triplets, pressureOffset, activeVsOffset, 1.0);
        addBlockEntriesToTripletsTranspose(JG_Matrix, triplets, pressureOffset, reducedVsOffset, 1.0);

        addBlockEntriesToTripletsTranspose(Dt_Matrix, triplets, stressOffset, activeVsOffset, 1.0);
        addBlockEntriesToTripletsTranspose(JDt_Matrix, triplets, stressOffset, reducedVsOffset, 1.0);
        addBlockEntriesToTriplets(uInv_Matrix, triplets, stressOffset, stressOffset, -0.5);

        addBlockEntriesToTripletsTranspose(JDtRed_Matrix, triplets, reducedStressOffset, reducedVsOffset, 1.0);
        addBlockEntriesToTriplets(uRed_Matrix, triplets, reducedStressOffset, reducedStressOffset, -0.5);

        A.setFromTriplets(triplets.begin(), triplets.end());
        //A.makeCompressed();
    }

    b = Vector::Zero(nSystemSize);
    {
        for (exint i = 0; i < nActiveVs; ++i)
            b(activeVsOffset + i) = invDt * activeRHSVector(i);
        for (exint i = 0; i < nReducedVs; ++i)
            b(reducedVsOffset + i) = invDt * reducedRHSVector(i);
        for (exint i = 0; i < nPressures; ++i)
            b(pressureOffset + i) = pressureRHSVector(i);
        for (exint i = 0; i < nStresses; ++i)
            b(stressOffset + i) = stressRHSVector(i);
    }

    guessVector = Vector::Zero(nSystemSize);
    {
        for (exint i = 0; i < nActiveVs; ++i)
            guessVector(activeVsOffset + i) = activeGuessVector(i);
        for (exint i = 0; i < nPressures; ++i)
            guessVector(pressureOffset + i) = pressureGuessVector(i);
        for (exint i = 0; i < nStresses; ++i)
            guessVector(stressOffset + i) = stressGuessVector(i);
        for (exint i = 0; i < nReducedVs; ++i)
            guessVector(reducedVsOffset + i) = reducedGuessVector(i);
    }

    solutionVector = Vector::Zero(nSystemSize);
}

void
HDK_PolyStokes::Solver::assembleSystemVelocityPressure()
{
    // create offsets
    exint activeVsOffset = 0;
    exint reducedVsOffset = nActiveVs;
    exint pressureOffset = nActiveVs + nReducedVs;
    nSystemSize = nActiveVs + nReducedVs + nPressures;

    // first assemble the reduced matrices and vectors
    //assembleReducedMassBlock();
    //assembleReducedViscosityBlock();
    assembleReducedCombinedBlock();
    assembleReducedRHSVector();
    assembleVMatrices();

    A.resize(nSystemSize, nSystemSize);
    // now assemble the entire thing together
    {
        Triplets triplets;
        exint sparseNonzeroElems =
            Mc_Matrix.nonZeros()
            + Mr_plus_2JDtuDJ_Matrix.nonZeros()
            + 2 * G_Matrix.nonZeros()
            + 2 * Dt_Matrix.nonZeros()
            + 2 * JG_Matrix.nonZeros()
            + 2 * JDt_Matrix.nonZeros()
            + u_Matrix.nonZeros();
        triplets.reserve(sparseNonzeroElems);

        // cwiseInverse is ok here since we know u is diagonal
        addBlockEntriesToTriplets(Mc_Matrix, triplets, activeVsOffset, activeVsOffset, invDt);
        addBlockEntriesToTriplets(V_Matrix, triplets, activeVsOffset, activeVsOffset, -1.);
        addBlockEntriesToTriplets(VJt_Matrix, triplets, activeVsOffset, reducedVsOffset, -1.);
        addBlockEntriesToTriplets(G_Matrix, triplets, activeVsOffset, pressureOffset, 1.);

        addBlockEntriesToTripletsTranspose(VJt_Matrix, triplets, reducedVsOffset, activeVsOffset, -1.);
        addBlockEntriesToTriplets(Mr_plus_2JDtuDJ_Matrix, triplets, reducedVsOffset, reducedVsOffset, 1.);          // includes only internal stresses
        addBlockEntriesToTriplets(JVJt_Matrix, triplets, reducedVsOffset, reducedVsOffset, -1.);                    // add in the boundary stresses
        addBlockEntriesToTriplets(JG_Matrix, triplets, reducedVsOffset, pressureOffset, 1.);

        addBlockEntriesToTripletsTranspose(G_Matrix, triplets, pressureOffset, activeVsOffset, 1.);
        addBlockEntriesToTripletsTranspose(JG_Matrix, triplets, pressureOffset, reducedVsOffset, 1.);

        A.setFromTriplets(triplets.begin(), triplets.end());
        //A.makeCompressed();
    }

    b = Vector::Zero(nSystemSize);
    {
        //Vector stressRHSContribution = Dt_Matrix * (2. * u_Matrix) * stressRHSVector;
        for (exint i = 0; i < nActiveVs; ++i)
            b(activeVsOffset + i) = invDt * activeRHSVector(i);// +stressRHSContribution(i);
        for (exint i = 0; i < nReducedVs; ++i)
            b(reducedVsOffset + i) = invDt * reducedRHSVector(i);
    }

    guessVector = Vector::Zero(nSystemSize);
    {
        for (exint i = 0; i < nActiveVs; ++i)
            guessVector(activeVsOffset + i) = activeGuessVector(i);
        for (exint i = 0; i < nReducedVs; ++i)
            guessVector(reducedVsOffset + i) = reducedGuessVector(i);
        for (exint i = 0; i < nPressures; ++i)
            guessVector(pressureOffset + i) = pressureGuessVector(i);
    }

    solutionVector = Vector::Zero(nSystemSize);
}

void
HDK_PolyStokes::Solver::assembleSystemPressureStress()
{
    // create offsets
    exint pressureOffset = 0;
    exint stressOffset = nPressures;
    nSystemSize = nPressures + nStresses;

    // first assemble the reduced matrices and vectors
    assembleReducedMassBlock();
    assembleReducedCombinedBlock();
    assembleReducedRHSVector();
    assembleReducedInvertedBlock();

    A.resize(nSystemSize, nSystemSize);
    // now assemble the entire thing together
    {
        Triplets triplets;
        // todo figure out a better estimate for this
        exint sparseNonzeroElems =
            McInv_Matrix.nonZeros()
            // + Mr_plus_2JDtuDJ_Matrix.nonZeros()
            + 2 * G_Matrix.nonZeros()
            + 2 * Dt_Matrix.nonZeros()
            // + 2 * JG_Matrix.nonZeros()
            // + 2 * JDt_Matrix.nonZeros()
            + uInv_Matrix.nonZeros();
        triplets.reserve(sparseNonzeroElems);

        // cwiseInverse is ok here since we know Mc is diagonal
        SparseMatrix A11_Matrix = -dt * G_Matrix.transpose() * McInv_Matrix * G_Matrix
            - JG_Matrix.transpose() * Inv_Mr_plus_2JDtuDJ_Matrix * JG_Matrix;
        SparseMatrix A12_Matrix = -dt * G_Matrix.transpose() * McInv_Matrix * Dt_Matrix
            - JG_Matrix.transpose() * Inv_Mr_plus_2JDtuDJ_Matrix * JDt_Matrix;
        SparseMatrix A21_Matrix = -dt * Dt_Matrix.transpose() * McInv_Matrix * G_Matrix
            - JDt_Matrix.transpose() * Inv_Mr_plus_2JDtuDJ_Matrix * JG_Matrix;
        SparseMatrix A22_Matrix = -dt * Dt_Matrix.transpose() * McInv_Matrix * Dt_Matrix
            - JDt_Matrix.transpose() * Inv_Mr_plus_2JDtuDJ_Matrix * JDt_Matrix
            - 0.5 * uInv_Matrix;

        addBlockEntriesToTriplets(A11_Matrix, triplets, pressureOffset, pressureOffset, 1.);
        addBlockEntriesToTriplets(A12_Matrix, triplets, pressureOffset, stressOffset, 1.);
        addBlockEntriesToTriplets(A21_Matrix, triplets, stressOffset, pressureOffset, 1.);
        //addBlockEntriesToTripletsTranspose(A12_Matrix, triplets, stressOffset, pressureOffset, 1.);
        addBlockEntriesToTriplets(A22_Matrix, triplets, stressOffset, stressOffset, 1.);

        A.setFromTriplets(triplets.begin(), triplets.end());
        //A.makeCompressed();
    }

    b = Vector::Zero(nSystemSize);
    {
        Vector b1_Vector = -G_Matrix.transpose() * McInv_Matrix * activeRHSVector
            - invDt * JG_Matrix.transpose() * Inv_Mr_plus_2JDtuDJ_Matrix * reducedRHSVector;
        Vector b2_Vector = -Dt_Matrix.transpose() * McInv_Matrix * activeRHSVector
            - invDt * JDt_Matrix.transpose() * Inv_Mr_plus_2JDtuDJ_Matrix * reducedRHSVector;

        /*
        for (exint i = 0; i < nPressures; ++i)
            b(pressureOffset + i) = b1_Vector(i);
        for (exint i = 0; i < nStresses; ++i)
            b(stressOffset + i) = b2_Vector(i);
        */

        for (exint i = 0; i < nPressures; ++i)
            b(pressureOffset + i) = b1_Vector(i) + pressureRHSVector(i);
        for (exint i = 0; i < nStresses; ++i)
            b(stressOffset + i) = b2_Vector(i) + stressRHSVector(i);
    }

    guessVector = Vector::Zero(nSystemSize);
    {
        for (exint i = 0; i < nPressures; ++i)
            guessVector(pressureOffset + i) = pressureGuessVector(i);
        for (exint i = 0; i < nStresses; ++i)
            guessVector(stressOffset + i) = stressGuessVector(i);
    }

    solutionVector = Vector::Zero(nSystemSize);
}

void
HDK_PolyStokes::Solver::assembleSystemPressureStressFactored()
{
    // create offsets
    exint pressureOffset = 0;
    exint stressOffset = nPressures;
    nSystemSize = nPressures + nStresses;

    // first assemble the reduced matrices and vectors
    assembleReducedMassBlock();
    assembleReducedCombinedBlock();
    assembleReducedRHSVector();
    assembleReducedInvertedBlock();

    A.resize(nSystemSize, nSystemSize);

    b = Vector::Zero(nSystemSize);
    {
        Vector b1_Vector = -G_Matrix.transpose() * McInv_Matrix * activeRHSVector
            - invDt * JG_Matrix.transpose() * Inv_Mr_plus_2JDtuDJ_Matrix * reducedRHSVector;
        Vector b2_Vector = -Dt_Matrix.transpose() * McInv_Matrix * activeRHSVector
            - invDt * JDt_Matrix.transpose() * Inv_Mr_plus_2JDtuDJ_Matrix * reducedRHSVector;

        for (exint i = 0; i < nPressures; ++i)
            b(pressureOffset + i) = b1_Vector(i) + pressureRHSVector(i);
        for (exint i = 0; i < nStresses; ++i)
            b(stressOffset + i) = b2_Vector(i) + stressRHSVector(i);
    }

    guessVector = Vector::Zero(nSystemSize);
    {
        for (exint i = 0; i < nPressures; ++i)
            guessVector(pressureOffset + i) = pressureGuessVector(i);
        for (exint i = 0; i < nStresses; ++i)
            guessVector(stressOffset + i) = stressGuessVector(i);
    }

    solutionVector = Vector::Zero(nSystemSize);
}