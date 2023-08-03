#include "HDK_PolyStokesSolver.h"

void
HDK_PolyStokes::Solver::assembleReducedMassBlock()
{
    Mr_Matrix.resize(nReducedVs, nReducedVs);
    Mr_Matrix.setZero();
    Mr_Matrix.data().squeeze();

    tbb::enumerable_thread_specific<Triplets> parallelMassMatrixElements;
    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            auto& localMassMatrixElements = parallelMassMatrixElements.local();
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
            {
                for (exint i = 0; i < REDUCED_DOF; ++i)
                    for (exint j = 0; j < REDUCED_DOF; ++j)
                        localMassMatrixElements.push_back(
                            Eigen::Triplet<SolveReal>(
                                REDUCED_DOF * interiorRegion + i,
                                REDUCED_DOF * interiorRegion + j,
                                reducedMassMatrices[interiorRegion](i, j))
                        );
            }
        });

    exint listSize = 0;

    parallelMassMatrixElements.combine_each([&](const Triplets& localMassMatrixElements)
        {
            listSize += localMassMatrixElements.size();
        });

    Triplets reducedMassMatrixElements;
    reducedMassMatrixElements.reserve(REDUCED_DOF * REDUCED_DOF * myInteriorRegionCount);
    reducedMassMatrixElements.reserve(listSize);

    parallelMassMatrixElements.combine_each([&](const Triplets& localMassMatrixElements)
        {
            reducedMassMatrixElements.insert(
                reducedMassMatrixElements.end(),
                localMassMatrixElements.begin(),
                localMassMatrixElements.end());
        });

    Mr_Matrix.setFromTriplets(reducedMassMatrixElements.begin(), reducedMassMatrixElements.end());
    //Mr_Matrix.makeCompressed();
}

void
HDK_PolyStokes::Solver::assembleReducedMassBlockInverse()
{
    MrInv_Matrix.resize(nReducedVs, nReducedVs);
    MrInv_Matrix.setZero();
    MrInv_Matrix.data().squeeze();

    tbb::enumerable_thread_specific<Triplets> parallelMassMatrixElements;
    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            auto& localMassMatrixElements = parallelMassMatrixElements.local();
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
            {
                ReducedMatrix localReducedMass = reducedMassMatrices[interiorRegion];
                ReducedMatrix invLocalReducedMass = localReducedMass.inverse();

                for (exint i = 0; i < REDUCED_DOF; ++i)
                    for (exint j = 0; j < REDUCED_DOF; ++j)
                        localMassMatrixElements.push_back(
                            Eigen::Triplet<SolveReal>(
                                REDUCED_DOF * interiorRegion + i,
                                REDUCED_DOF * interiorRegion + j,
                                invLocalReducedMass(i, j))
                        );
            }
        });

    exint listSize = 0;

    parallelMassMatrixElements.combine_each([&](const Triplets& localMassMatrixElements)
        {
            listSize += localMassMatrixElements.size();
        });

    Triplets reducedMassMatrixElements;
    reducedMassMatrixElements.reserve(REDUCED_DOF * REDUCED_DOF * myInteriorRegionCount);
    reducedMassMatrixElements.reserve(listSize);

    parallelMassMatrixElements.combine_each([&](const Triplets& localMassMatrixElements)
        {
            reducedMassMatrixElements.insert(
                reducedMassMatrixElements.end(),
                localMassMatrixElements.begin(),
                localMassMatrixElements.end());
        });

    MrInv_Matrix.setFromTriplets(reducedMassMatrixElements.begin(), reducedMassMatrixElements.end());
    //MrInv_Matrix.makeCompressed();
}

void
HDK_PolyStokes::Solver::assembleReducedViscosityBlock()
{
    JDtuDJ_Matrix.resize(nReducedVs, nReducedVs);
    JDtuDJ_Matrix.setZero();
    JDtuDJ_Matrix.data().squeeze();

    tbb::enumerable_thread_specific<Triplets> parallelViscosityMatrixElements;
    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            auto& localViscosityMatrixElements = parallelViscosityMatrixElements.local();
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
            {
                for (exint i = 0; i < REDUCED_DOF; ++i)
                    for (exint j = 0; j < REDUCED_DOF; ++j)
                        localViscosityMatrixElements.push_back(
                            Eigen::Triplet<SolveReal>(
                                REDUCED_DOF * interiorRegion + i,
                                REDUCED_DOF * interiorRegion + j,
                                reducedViscosityMatrices[interiorRegion](i, j))
                        );
            }
        });

    exint listSize = 0;

    parallelViscosityMatrixElements.combine_each([&](const Triplets& localViscosityMatrixElements)
        {
            listSize += localViscosityMatrixElements.size();
        });

    Triplets reducedViscosityMatrixElements;
    reducedViscosityMatrixElements.reserve(REDUCED_DOF * REDUCED_DOF * myInteriorRegionCount);
    reducedViscosityMatrixElements.reserve(listSize);

    parallelViscosityMatrixElements.combine_each([&](const Triplets& localViscosityMatrixElements)
        {
            reducedViscosityMatrixElements.insert(
                reducedViscosityMatrixElements.end(),
                localViscosityMatrixElements.begin(),
                localViscosityMatrixElements.end());
        });

    JDtuDJ_Matrix.setFromTriplets(reducedViscosityMatrixElements.begin(), reducedViscosityMatrixElements.end());
    //JDtuDJ_Matrix.makeCompressed();
}

void
HDK_PolyStokes::Solver::assembleReducedCombinedBlock()
{
    Mr_plus_2JDtuDJ_Matrix.resize(nReducedVs, nReducedVs);
    Mr_plus_2JDtuDJ_Matrix.setZero();
    Mr_plus_2JDtuDJ_Matrix.data().squeeze();

    tbb::enumerable_thread_specific<Triplets> parallelReducedMatrixElements;
    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            auto& localReducedMatrixElements = parallelReducedMatrixElements.local();
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
            {
                for (exint i = 0; i < REDUCED_DOF; ++i)
                    for (exint j = 0; j < REDUCED_DOF; ++j)
                        localReducedMatrixElements.push_back(
                            Eigen::Triplet<SolveReal>(
                                REDUCED_DOF * interiorRegion + i,
                                REDUCED_DOF * interiorRegion + j,
                                invDt * reducedMassMatrices[interiorRegion](i, j) + 2. * reducedViscosityMatrices[interiorRegion](i, j)
                                )
                        );
            }
        });

    exint listSize = 0;

    parallelReducedMatrixElements.combine_each([&](const Triplets& localReducedMatrixElements)
        {
            listSize += localReducedMatrixElements.size();
        });

    Triplets reducedMatrixElements;
    reducedMatrixElements.reserve(REDUCED_DOF * REDUCED_DOF * myInteriorRegionCount);
    reducedMatrixElements.reserve(listSize);

    parallelReducedMatrixElements.combine_each([&](const Triplets& localReducedMatrixElements)
        {
            reducedMatrixElements.insert(
                reducedMatrixElements.end(),
                localReducedMatrixElements.begin(),
                localReducedMatrixElements.end());
        });

    Mr_plus_2JDtuDJ_Matrix.setFromTriplets(reducedMatrixElements.begin(), reducedMatrixElements.end());
    //Mr_plus_2JDtuDJ_Matrix.makeCompressed();
}

void
HDK_PolyStokes::Solver::assembleReducedInvertedBlock()
{
    Inv_Mr_plus_2JDtuDJ_Matrix.resize(nReducedVs, nReducedVs);
    Inv_Mr_plus_2JDtuDJ_Matrix.setZero();
    Inv_Mr_plus_2JDtuDJ_Matrix.data().squeeze();

    tbb::enumerable_thread_specific<Triplets> parallelReducedMatrixElements;
    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            auto& localReducedMatrixElements = parallelReducedMatrixElements.local();
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
            {
                ReducedMatrix localReducedMatrix = invDt * reducedMassMatrices[interiorRegion] + 2. * reducedViscosityMatrices[interiorRegion];
                ReducedMatrix invLocalReducedMatrix = localReducedMatrix.inverse();

                for (exint i = 0; i < REDUCED_DOF; ++i)
                    for (exint j = 0; j < REDUCED_DOF; ++j)
                        localReducedMatrixElements.push_back(
                            Eigen::Triplet<SolveReal>(
                                REDUCED_DOF * interiorRegion + i,
                                REDUCED_DOF * interiorRegion + j,
                                invLocalReducedMatrix(i, j)
                                )
                        );
            }
        });

    exint listSize = 0;

    parallelReducedMatrixElements.combine_each([&](const Triplets& localReducedMatrixElements)
        {
            listSize += localReducedMatrixElements.size();
        });

    Triplets reducedMatrixElements;
    reducedMatrixElements.reserve(REDUCED_DOF * REDUCED_DOF * myInteriorRegionCount);
    reducedMatrixElements.reserve(listSize);

    parallelReducedMatrixElements.combine_each([&](const Triplets& localReducedMatrixElements)
        {
            reducedMatrixElements.insert(
                reducedMatrixElements.end(),
                localReducedMatrixElements.begin(),
                localReducedMatrixElements.end());
        });

    Inv_Mr_plus_2JDtuDJ_Matrix.setFromTriplets(reducedMatrixElements.begin(), reducedMatrixElements.end());
    //Inv_Mr_plus_2JDtuDJ_Matrix.makeCompressed();
}

void
HDK_PolyStokes::Solver::assembleReducedInvertedBlockIncludingBoundaryStresses()
{
    // note requires JVJt to exist (ie call assembleVMatrices before this)

    BInv_Matrix.resize(nReducedVs, nReducedVs);
    BInv_Matrix.setZero();
    BInv_Matrix.data().squeeze();

    tbb::enumerable_thread_specific<Triplets> parallelReducedMatrixElements;
    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            auto& localReducedMatrixElements = parallelReducedMatrixElements.local();
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
            {
                ReducedMatrix localReducedMatrix = invDt * reducedMassMatrices[interiorRegion]
                    + 2. * reducedViscosityMatrices[interiorRegion]
                    - JVJt_Matrix.block(REDUCED_DOF * interiorRegion, REDUCED_DOF * interiorRegion, REDUCED_DOF, REDUCED_DOF);
                ReducedMatrix invLocalReducedMatrix = localReducedMatrix.inverse();

                for (exint i = 0; i < REDUCED_DOF; ++i)
                    for (exint j = 0; j < REDUCED_DOF; ++j)
                        localReducedMatrixElements.push_back(
                            Eigen::Triplet<SolveReal>(
                                REDUCED_DOF * interiorRegion + i,
                                REDUCED_DOF * interiorRegion + j,
                                invLocalReducedMatrix(i, j)
                                )
                        );
            }
        });

    exint listSize = 0;

    parallelReducedMatrixElements.combine_each([&](const Triplets& localReducedMatrixElements)
        {
            listSize += localReducedMatrixElements.size();
        });

    Triplets reducedMatrixElements;
    reducedMatrixElements.reserve(REDUCED_DOF * REDUCED_DOF * myInteriorRegionCount);
    reducedMatrixElements.reserve(listSize);

    parallelReducedMatrixElements.combine_each([&](const Triplets& localReducedMatrixElements)
        {
            reducedMatrixElements.insert(
                reducedMatrixElements.end(),
                localReducedMatrixElements.begin(),
                localReducedMatrixElements.end());
        });

    BInv_Matrix.setFromTriplets(reducedMatrixElements.begin(), reducedMatrixElements.end());
    //BInv_Matrix.makeCompressed();
}

void
HDK_PolyStokes::Solver::assembleReducedUninvertedBlockIncludingBoundaryStresses()
{
    // note requires JVJt to exist (ie call assembleVMatrices before this)

    B_Matrix.resize(nReducedVs, nReducedVs);
    B_Matrix.setZero();
    B_Matrix.data().squeeze();

    tbb::enumerable_thread_specific<Triplets> parallelReducedMatrixElements;
    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            auto& localReducedMatrixElements = parallelReducedMatrixElements.local();
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
            {
                ReducedMatrix localReducedMatrix = invDt * reducedMassMatrices[interiorRegion]
                    + 2. * reducedViscosityMatrices[interiorRegion]
                    - JVJt_Matrix.block(REDUCED_DOF * interiorRegion, REDUCED_DOF * interiorRegion, REDUCED_DOF, REDUCED_DOF);
                ReducedMatrix invLocalReducedMatrix = localReducedMatrix;

                for (exint i = 0; i < REDUCED_DOF; ++i)
                    for (exint j = 0; j < REDUCED_DOF; ++j)
                        localReducedMatrixElements.push_back(
                            Eigen::Triplet<SolveReal>(
                                REDUCED_DOF * interiorRegion + i,
                                REDUCED_DOF * interiorRegion + j,
                                invLocalReducedMatrix(i, j)
                                )
                        );
            }
        });

    exint listSize = 0;

    parallelReducedMatrixElements.combine_each([&](const Triplets& localReducedMatrixElements)
        {
            listSize += localReducedMatrixElements.size();
        });

    Triplets reducedMatrixElements;
    reducedMatrixElements.reserve(REDUCED_DOF * REDUCED_DOF * myInteriorRegionCount);
    reducedMatrixElements.reserve(listSize);

    parallelReducedMatrixElements.combine_each([&](const Triplets& localReducedMatrixElements)
        {
            reducedMatrixElements.insert(
                reducedMatrixElements.end(),
                localReducedMatrixElements.begin(),
                localReducedMatrixElements.end());
        });

    B_Matrix.setFromTriplets(reducedMatrixElements.begin(), reducedMatrixElements.end());
    //B_Matrix.makeCompressed();
}

void
HDK_PolyStokes::Solver::assembleReducedRHSVector()
{
    reducedRHSVector = Vector::Zero(nReducedVs);

    for (int i = 0; i < myInteriorRegionCount; ++i)
    {
        ColumnVector localRHSVector = reducedMassMatrices[i] * reducedRegionBestFitVectors[i];
        for (int j = 0; j < REDUCED_DOF; ++j)
            reducedRHSVector(REDUCED_DOF * i + j) = localRHSVector(j);
    }
}

void
HDK_PolyStokes::Solver::assembleVMatrices()
{
    V_Matrix = Dt_Matrix * (-2. * u_Matrix) * Dt_Matrix.transpose();
    VJt_Matrix = Dt_Matrix * (-2. * u_Matrix) * JDt_Matrix.transpose();
    JVJt_Matrix = JDt_Matrix * (-2. * u_Matrix) * JDt_Matrix.transpose();
}
