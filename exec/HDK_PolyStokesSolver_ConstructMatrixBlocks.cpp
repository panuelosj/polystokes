#include "HDK_PolyStokesSolver.h"

// gas solver tools
#include <SIM/SIM_Object.h>
// multithreading
#include <UT/UT_ParallelUtil.h>
#include <tbb/tbb.h>

void
HDK_PolyStokes::Solver::constructMatrixBlocks()
{
    // first figure out all our dimensions
    nActiveVx = nFaceX; nActiveVy = nFaceY; nActiveVz = nFaceZ;
    nActiveVs = nActiveVx + nActiveVy + nActiveVz;
    nReducedVs = myInteriorRegionCount * REDUCED_DOF;
    nPressures = nCenter;
    nTxx = nCenter; nTyy = nCenter; nTzz = nCenter;
    nTxy = nEdgeXY; nTxz = nEdgeXZ; nTyz = nEdgeYZ;
    nStresses = nTxx + nTyy + nTzz + nTxy + nTxz + nTyz;
    nReducedStresses = myInteriorRegionCount * 6; // one for each T variable
    nTotalDOFs = nActiveVs + nReducedVs + nPressures + nStresses;

    // resize blocks to their appropriate sizes
    Mc_Matrix.resize(nActiveVs, nActiveVs);
    McInv_Matrix.resize(nActiveVs, nActiveVs);
    uInv_Matrix.resize(nStresses, nStresses);
    u_Matrix.resize(nStresses, nStresses);
    uInvRed_Matrix.resize(nReducedStresses, nReducedStresses);
    uRed_Matrix.resize(nReducedStresses, nReducedStresses);
    // off diagonal blocks
    G_Matrix.resize(nActiveVs, nPressures);
    Dt_Matrix.resize(nActiveVs, nStresses);
    JG_Matrix.resize(nReducedVs, nPressures);
    JDt_Matrix.resize(nReducedVs, nStresses);
    JDtRed_Matrix.resize(nReducedVs, nReducedStresses);
    // rhs vectors
    activeRHSVectorSparse.resize(nActiveVs, 1);
    pressureRHSVectorSparse.resize(nPressures, 1);
    stressRHSVectorSparse.resize(nStresses, 1);
    // save old velocities to use for guess vectors
    oldActiveVs.resize(nActiveVs, 1);


    // parallel list of triplets
    std::vector<Triplets> parallelActiveMassMatrixElements(myThreadCount);
    std::vector<Triplets> parallelActiveMassInverseMatrixElements(myThreadCount);
    std::vector<Triplets> parallelActiveStressInverseMatrixElements(myThreadCount);
    std::vector<Triplets> parallelActiveStressMatrixElements(myThreadCount);
    std::vector<Triplets> parallelReducedStressInverseMatrixElements(myThreadCount);
    std::vector<Triplets> parallelReducedStressMatrixElements(myThreadCount);
    std::vector<Triplets> parallelActiveGradientMatrixElements(myThreadCount);
    std::vector<Triplets> parallelReducedGradientMatrixElements(myThreadCount);
    std::vector<Triplets> parallelActiveDivergenceMatrixElements(myThreadCount);
    std::vector<Triplets> parallelReducedDivergenceMatrixElements(myThreadCount);
    std::vector<Triplets> parallelReducedInternalDivergenceMatrixElements(myThreadCount);
    std::vector<Triplets> parallelActiveRHSElements(myThreadCount);
    std::vector<Triplets> parallelPressureRHSElements(myThreadCount);
    std::vector<Triplets> parallelStressRHSElements(myThreadCount);
    std::vector<Triplets> parallelActiveOldVsElements(myThreadCount);
    buildMatrixBlocksByTriplets(
        parallelActiveMassMatrixElements,
        parallelActiveMassInverseMatrixElements,
        parallelActiveStressInverseMatrixElements,
        parallelActiveStressMatrixElements,
        parallelReducedStressInverseMatrixElements,
        parallelReducedStressMatrixElements,
        parallelActiveGradientMatrixElements,
        parallelReducedGradientMatrixElements,
        parallelActiveDivergenceMatrixElements,
        parallelReducedDivergenceMatrixElements,
        parallelReducedInternalDivergenceMatrixElements,
        parallelActiveRHSElements,
        parallelPressureRHSElements,
        parallelStressRHSElements,
        parallelActiveOldVsElements
    );

    // now compile
    // Compile RHS vector
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelActiveRHSElements[thread].size();

        Triplets activeRHSElements;
        activeRHSElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            activeRHSElements.insert(activeRHSElements.end(), parallelActiveRHSElements[thread].begin(), parallelActiveRHSElements[thread].end());

        activeRHSVectorSparse.setFromTriplets(activeRHSElements.begin(), activeRHSElements.end());
        activeRHSVector = Vector(activeRHSVectorSparse);
    }
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelPressureRHSElements[thread].size();

        Triplets pressureRHSElements;
        pressureRHSElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            pressureRHSElements.insert(pressureRHSElements.end(), parallelPressureRHSElements[thread].begin(), parallelPressureRHSElements[thread].end());

        pressureRHSVectorSparse.setFromTriplets(pressureRHSElements.begin(), pressureRHSElements.end());
        pressureRHSVector = Vector(pressureRHSVectorSparse);
    }
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelStressRHSElements[thread].size();

        Triplets stressRHSElements;
        stressRHSElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            stressRHSElements.insert(stressRHSElements.end(), parallelStressRHSElements[thread].begin(), parallelStressRHSElements[thread].end());

        stressRHSVectorSparse.setFromTriplets(stressRHSElements.begin(), stressRHSElements.end());
        stressRHSVector = Vector(stressRHSVectorSparse);
    }

    // Compile old velocity vector
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelActiveOldVsElements[thread].size();

        Triplets activeOldVsElements;
        activeOldVsElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            activeOldVsElements.insert(activeOldVsElements.end(), parallelActiveOldVsElements[thread].begin(), parallelActiveOldVsElements[thread].end());

        oldActiveVs.setFromTriplets(activeOldVsElements.begin(), activeOldVsElements.end());
        //oldActiveVs.makeCompressed();
    }
    // Compile mass matrix
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelActiveMassMatrixElements[thread].size();

        Triplets matrixElements;
        matrixElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            matrixElements.insert(matrixElements.end(), parallelActiveMassMatrixElements[thread].begin(), parallelActiveMassMatrixElements[thread].end());

        Mc_Matrix.setFromTriplets(matrixElements.begin(), matrixElements.end());
        //totalMatrix.makeCompressed();
    }
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelActiveMassInverseMatrixElements[thread].size();

        Triplets matrixElements;
        matrixElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            matrixElements.insert(matrixElements.end(), parallelActiveMassInverseMatrixElements[thread].begin(), parallelActiveMassInverseMatrixElements[thread].end());

        McInv_Matrix.setFromTriplets(matrixElements.begin(), matrixElements.end());
        //totalMatrix.makeCompressed();
    }
    // Compile stress diagonal matrix
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelActiveStressInverseMatrixElements[thread].size();

        Triplets matrixElements;
        matrixElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            matrixElements.insert(matrixElements.end(), parallelActiveStressInverseMatrixElements[thread].begin(), parallelActiveStressInverseMatrixElements[thread].end());

        uInv_Matrix.setFromTriplets(matrixElements.begin(), matrixElements.end());
        //totalMatrix.makeCompressed();
    }
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelActiveStressMatrixElements[thread].size();

        Triplets matrixElements;
        matrixElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            matrixElements.insert(matrixElements.end(), parallelActiveStressMatrixElements[thread].begin(), parallelActiveStressMatrixElements[thread].end());

        u_Matrix.setFromTriplets(matrixElements.begin(), matrixElements.end());
        //totalMatrix.makeCompressed();
    }
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelReducedStressInverseMatrixElements[thread].size();

        Triplets matrixElements;
        matrixElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            matrixElements.insert(matrixElements.end(), parallelReducedStressInverseMatrixElements[thread].begin(), parallelReducedStressInverseMatrixElements[thread].end());

        uInvRed_Matrix.setFromTriplets(matrixElements.begin(), matrixElements.end());
    }
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelReducedStressMatrixElements[thread].size();

        Triplets matrixElements;
        matrixElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            matrixElements.insert(matrixElements.end(), parallelReducedStressMatrixElements[thread].begin(), parallelReducedStressMatrixElements[thread].end());

        uRed_Matrix.setFromTriplets(matrixElements.begin(), matrixElements.end());
    }
    // Compile active gradient matrix
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelActiveGradientMatrixElements[thread].size();

        Triplets matrixElements;
        matrixElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            matrixElements.insert(matrixElements.end(), parallelActiveGradientMatrixElements[thread].begin(), parallelActiveGradientMatrixElements[thread].end());

        G_Matrix.setFromTriplets(matrixElements.begin(), matrixElements.end());
    }
    // Compile reduced gradient matrix
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelReducedGradientMatrixElements[thread].size();

        Triplets matrixElements;
        matrixElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            matrixElements.insert(matrixElements.end(), parallelReducedGradientMatrixElements[thread].begin(), parallelReducedGradientMatrixElements[thread].end());

        JG_Matrix.setFromTriplets(matrixElements.begin(), matrixElements.end());
    }
    // Compile active divergence matrix
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelActiveDivergenceMatrixElements[thread].size();

        Triplets matrixElements;
        matrixElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            matrixElements.insert(matrixElements.end(), parallelActiveDivergenceMatrixElements[thread].begin(), parallelActiveDivergenceMatrixElements[thread].end());

        Dt_Matrix.setFromTriplets(matrixElements.begin(), matrixElements.end());
        //totalMatrix.makeCompressed();
    }
    // Compile reduced divergence matrix
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelReducedDivergenceMatrixElements[thread].size();

        Triplets matrixElements;
        matrixElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            matrixElements.insert(matrixElements.end(), parallelReducedDivergenceMatrixElements[thread].begin(), parallelReducedDivergenceMatrixElements[thread].end());

        JDt_Matrix.setFromTriplets(matrixElements.begin(), matrixElements.end());
    }
    {
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelReducedInternalDivergenceMatrixElements[thread].size();

        Triplets matrixElements;
        matrixElements.reserve(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
            matrixElements.insert(matrixElements.end(), parallelReducedInternalDivergenceMatrixElements[thread].begin(), parallelReducedInternalDivergenceMatrixElements[thread].end());

        JDtRed_Matrix.setFromTriplets(matrixElements.begin(), matrixElements.end());
    }
}

void
HDK_PolyStokes::Solver::buildMatrixBlocksByTriplets(
    std::vector<Triplets>& parallelActiveMassMatrixElements,
    std::vector<Triplets>& parallelActiveMassInverseMatrixElements,
    std::vector<Triplets>& parallelActiveStressInverseMatrixElements,
    std::vector<Triplets>& parallelActiveStressMatrixElements,
    std::vector<Triplets>& parallelReducedStressInverseMatrixElements,
    std::vector<Triplets>& parallelReducedStressMatrixElements,
    std::vector<Triplets>& parallelActiveGradientMatrixElements,
    std::vector<Triplets>& parallelReducedGradientMatrixElements,
    std::vector<Triplets>& parallelActiveDivergenceMatrixElements,
    std::vector<Triplets>& parallelReducedDivergenceMatrixElements,
    std::vector<Triplets>& parallelReducedInternalDivergenceMatrixElements,
    std::vector<Triplets>& parallelActiveRHSElements,
    std::vector<Triplets>& parallelPressureRHSElements,
    std::vector<Triplets>& parallelStressRHSElements,
    std::vector<Triplets>& parallelActiveOldVsElements
)
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::faceToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::faceToEdgeMap;
    using SIM::FieldUtils::edgeToFaceMap;

    // setup stencils for mass, G, Dt
    for (int faceAxis : {0, 1, 2})
    {
        UT_Interrupt* boss = UTgetInterrupt();

        UT_ThreadedAlgorithm buildLinearSystemAlgorithm;
        buildLinearSystemAlgorithm.run([&](const UT_JobInfo& info)
            {
                Triplets& localActiveMassMatrixElements = parallelActiveMassMatrixElements[info.job()];
                Triplets& localActiveMassInverseMatrixElements = parallelActiveMassInverseMatrixElements[info.job()];
                Triplets& localActiveGradientMatrixElements = parallelActiveGradientMatrixElements[info.job()];
                Triplets& localReducedGradientMatrixElements = parallelReducedGradientMatrixElements[info.job()];
                Triplets& localActiveDivergenceMatrixElements = parallelActiveDivergenceMatrixElements[info.job()];
                Triplets& localReducedDivergenceMatrixElements = parallelReducedDivergenceMatrixElements[info.job()];
                Triplets& localReducedInternalDivergenceMatrixElements = parallelReducedInternalDivergenceMatrixElements[info.job()];
                Triplets& localActiveRHSElements = parallelActiveRHSElements[info.job()];
                Triplets& localPressureRHSElements = parallelPressureRHSElements[info.job()];
                Triplets& localStressRHSElements = parallelStressRHSElements[info.job()];
                Triplets& localActiveOldVsElements = parallelActiveOldVsElements[info.job()];

                UT_VoxelArrayIteratorI vit;
                vit.setConstArray(faceLabels(faceAxis)->field());
                vit.splitByTile(info);

                UT_VoxelTileIteratorI vitt;
                for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
                {
                    if (boss->opInterrupt())
                        break;

                    if (!vit.isTileConstant() || isSolved(vit.getValue()))
                    {
                        vitt.setTile(vit);
                        for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                        {
                            //int j = faceAxis;
                            UT_Vector3I face(vitt.x(), vitt.y(), vitt.z());

                            MaterialLabels selfLabel = (MaterialLabels)getFieldValue(*faceLabels(faceAxis), face);
                            exint selfActiveIndex = faceVelocityDOF(getFieldValue(*faceActiveIndices(faceAxis), face), faceAxis);
                            exint selfReducedIndex = getFieldValue(*faceReducedIndices(faceAxis), face);

                            SolveReal localDensity = getLocalDensity(face);
                            SolveReal volume =
                                (SolveReal)getFieldValue(*faceFluidWeights(faceAxis), face)
                                * (SolveReal)getFieldValue(*faceLiquidWeights(faceAxis), face);
                            volume = SYSclamp(volume, MINWEIGHT * MINWEIGHT, 1.0);
                            SolveReal localVelocity = getFieldValue(*myVelocityField->getField(faceAxis), face);

                            // setup diagonal mass matrix, rhs, and warm start guess
                            if (isActive(selfLabel))
                            {
                                localActiveMassMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    selfActiveIndex,
                                    selfActiveIndex,
                                    volume * localDensity)
                                );
                                localActiveMassInverseMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    selfActiveIndex,
                                    selfActiveIndex,
                                    1. / (volume * localDensity))
                                );
                                localActiveRHSElements.push_back(Eigen::Triplet<SolveReal>(
                                    selfActiveIndex,
                                    0,
                                    localVelocity * volume * localDensity)
                                );
                                localActiveOldVsElements.push_back(Eigen::Triplet<SolveReal>(
                                    selfActiveIndex,
                                    0,
                                    localVelocity)
                                );
                            }

                            // setup pressure stencils
                            if (isActive(selfLabel) || isReduced(selfLabel))
                            {
                                for (int gradientDirection : {0, 1})
                                {
                                    SolveReal gradientSign = (gradientDirection == 0) ? -1 : 1;
                                    UT_Vector3I cell = faceToCellMap(face, faceAxis, gradientDirection);
                                    // oob check
                                    if (cell[faceAxis] < 0 || cell[faceAxis] >= faceActiveIndices(faceAxis)->getVoxelRes()[faceAxis])
                                        continue;

                                    exint cellPressureIndex = getFieldValue(centerActiveIndices, cell);

                                    if (cellPressureIndex >= 0)
                                    {
                                        SolveReal coeff =
                                            (SolveReal)getFieldValue(*faceFluidWeights(faceAxis), face)
                                            * (SolveReal)getFieldValue(centerLiquidWeights, cell)
                                            * (SolveReal)invDx;
                                        SolveReal contribution = gradientSign * coeff;

                                        if (coeff <= 0.) continue;
                                        if (coeff > 0.)
                                        {
                                            if (isActive(selfLabel))
                                            {
                                                // pressureGradient
                                                localActiveGradientMatrixElements.push_back(Eigen::Triplet<SolveReal>(selfActiveIndex, cellPressureIndex, contribution));

                                                // todo fix the coeff constants in the solid bdd
                                                // solid boundary
                                                if (getFieldValue(centerFluidWeights, cell) < 1.)
                                                {
                                                    SolveReal solidCenterVolume = 1. - getFieldValue(centerFluidWeights, cell);
                                                    SolveReal liquidCenterVolume = getFieldValue(centerLiquidWeights, cell);
                                                    SolveReal solidCoeff = liquidCenterVolume * solidCenterVolume * invDx;
                                                    SolveReal solidContribution = gradientSign * coeff;
                                                    SolveReal svel = getFieldValue(*myCollisionVelocityField->getField(faceAxis), face);
                                                    localPressureRHSElements.push_back(Eigen::Triplet<SolveReal>(cellPressureIndex, 0, -1. * solidContribution * svel));
                                                }
                                                if (getFieldValue(*faceFluidWeights(faceAxis), face) < 1.)
                                                { // RHS Contribution
                                                    SolveReal solidFaceVolume = 1. - getFieldValue(*faceFluidWeights(faceAxis), face);
                                                    SolveReal liquidCenterVolume = getFieldValue(centerLiquidWeights, cell);
                                                    SolveReal solidCoeff = liquidCenterVolume * solidFaceVolume * invDx;
                                                    SolveReal solidContribution = gradientSign * coeff;
                                                    SolveReal svel = getFieldValue(*myCollisionVelocityField->getField(faceAxis), face);
                                                    localPressureRHSElements.push_back(Eigen::Triplet<SolveReal>(cellPressureIndex, 0, solidContribution * svel));
                                                }
                                            }
                                            else if (isReduced(selfLabel))
                                            {
                                                // setup offsets
                                                UT_Vector3T<SolveReal> selfOffset;
                                                UT_Vector3T<SolveReal> facePos(face);
                                                facePos[faceAxis] -= 0.5;
                                                facePos *= dx;
                                                selfOffset = facePos - reducedRegionCOM[selfReducedIndex];

                                                ColumnVector colVec = buildConversionCoefficients(selfOffset, faceAxis);

                                                for (int n = 0; n < REDUCED_DOF; n++)
                                                    localReducedGradientMatrixElements.push_back(Eigen::Triplet<SolveReal>(REDUCED_DOF * selfReducedIndex + n, cellPressureIndex, contribution * colVec(n)));
                                            }
                                        }
                                    }
                                }
                            }

                            // setup stress stencils
                            if (isActive(selfLabel) || isReduced(selfLabel))
                            {
                                // centers
                                for (int divergenceDirection : {0, 1})
                                {
                                    SolveReal divergenceSign = (divergenceDirection == 0) ? -1 : 1;
                                    UT_Vector3I cell = faceToCellMap(face, faceAxis, divergenceDirection);
                                    // oob check
                                    if (cell[faceAxis] < 0 || cell[faceAxis] >= faceActiveIndices(faceAxis)->getVoxelRes()[faceAxis])
                                        continue;

                                    // we need to make sure we translate the axis DOF into the index of the entire stress DOFs
                                    MaterialLabels cellLabel = (MaterialLabels)getFieldValue(centerLabels, cell);
                                    exint cellStressIndex = centerStressDOF(getFieldValue(centerActiveIndices, cell), faceAxis);

                                    if (isActive(cellLabel))
                                    {
                                        SolveReal coeff =
                                            (SolveReal)getFieldValue(*faceFluidWeights(faceAxis), face)
                                            * (SolveReal)getFieldValue(centerLiquidWeights, cell)
                                            * (SolveReal)invDx;
                                        SolveReal contribution = -1. * divergenceSign * coeff;

                                        if (coeff <= 0.) continue;
                                        if (coeff > 0.)
                                        {
                                            if (isActive(selfLabel))
                                            {
                                                localActiveDivergenceMatrixElements.push_back(Eigen::Triplet<SolveReal>(selfActiveIndex, cellStressIndex, contribution));

                                                // solid boundary
                                                if (getFieldValue(centerFluidWeights, cell) < 1.)
                                                {
                                                    SolveReal solidCenterVolume = 1. - getFieldValue(centerFluidWeights, cell);
                                                    SolveReal liquidCenterVolume = getFieldValue(centerLiquidWeights, cell);
                                                    SolveReal solidCoeff = liquidCenterVolume * solidCenterVolume * invDx;
                                                    SolveReal solidContribution = divergenceSign * coeff;
                                                    SolveReal svel = getFieldValue(*myCollisionVelocityField->getField(faceAxis), face);
                                                    localStressRHSElements.push_back(Eigen::Triplet<SolveReal>(cellStressIndex, 0, -1. * solidContribution * svel));
                                                }
                                                if (getFieldValue(*faceFluidWeights(faceAxis), face) < 1.)
                                                { // RHS Contribution
                                                    SolveReal solidFaceVolume = 1. - getFieldValue(*faceFluidWeights(faceAxis), face);
                                                    SolveReal liquidCenterVolume = getFieldValue(centerLiquidWeights, cell);
                                                    SolveReal solidCoeff = liquidCenterVolume * solidFaceVolume * invDx;
                                                    SolveReal solidContribution = divergenceSign * coeff;
                                                    SolveReal svel = getFieldValue(*myCollisionVelocityField->getField(faceAxis), face);
                                                    localStressRHSElements.push_back(Eigen::Triplet<SolveReal>(cellStressIndex, 0, solidContribution * svel));
                                                }
                                            }
                                            else if (isReduced(selfLabel))
                                            {
                                                // setup offsets
                                                UT_Vector3T<SolveReal> selfOffset;
                                                UT_Vector3T<SolveReal> facePos(face);
                                                facePos[faceAxis] -= 0.5;
                                                facePos *= dx;
                                                selfOffset = facePos - reducedRegionCOM[selfReducedIndex];

                                                ColumnVector colVec = buildConversionCoefficients(selfOffset, faceAxis);

                                                for (int n = 0; n < REDUCED_DOF; n++)
                                                    localReducedDivergenceMatrixElements.push_back(Eigen::Triplet<SolveReal>(REDUCED_DOF * selfReducedIndex + n, cellStressIndex, contribution * colVec(n)));
                                            }
                                        }
                                    }
                                    if (isReduced(cellLabel))
                                    {
                                        exint reducedCellIndex = getFieldValue(centerReducedIndices, cell);
                                        exint reducedCellStressIndex = reducedCenterStressDOF(reducedCellIndex, faceAxis);
                                        SolveReal contribution = -1. * invDx * divergenceSign;

                                        if (isReduced(selfLabel))
                                        {
                                            // setup offsets
                                            UT_Vector3T<SolveReal> selfOffset;
                                            UT_Vector3T<SolveReal> facePos(face);
                                            facePos[faceAxis] -= 0.5;
                                            facePos *= dx;
                                            selfOffset = facePos - reducedRegionCOM[selfReducedIndex];

                                            ColumnVector colVec = buildConversionCoefficients(selfOffset, faceAxis);

                                            for (int n = 0; n < REDUCED_DOF; n++)
                                                localReducedInternalDivergenceMatrixElements.push_back(Eigen::Triplet<SolveReal>(REDUCED_DOF * selfReducedIndex + n, reducedCellStressIndex, contribution * colVec(n)));
                                        }
                                    }
                                }

                                // edges
                                for (int edgeAxis : {0, 1, 2})
                                {
                                    if (edgeAxis == faceAxis) continue;

                                    for (int divergenceDirection : {0, 1})
                                    {
                                        SolveReal divergenceSign = (divergenceDirection == 0) ? -1 : 1;
                                        UT_Vector3I edge = faceToEdgeMap(face, faceAxis, edgeAxis, divergenceDirection);

                                        // we need to make sure we translate the axis DOF into the index of the entire stress DOFs
                                        MaterialLabels edgeLabel = (MaterialLabels)getFieldValue(*edgeLabels(edgeAxis), edge);
                                        exint edgeStressIndex = edgeStressDOF(getFieldValue(*edgeActiveIndices(edgeAxis), edge), edgeAxis);

                                        if (isActive(edgeLabel))
                                        {
                                            SolveReal coeff =
                                                (SolveReal)getFieldValue(*faceFluidWeights(faceAxis), face)
                                                * (SolveReal)getFieldValue(*edgeLiquidWeights(edgeAxis), edge)
                                                * (SolveReal)invDx;
                                            SolveReal contribution = -1. * divergenceSign * coeff;

                                            if (coeff <= 0.) continue;
                                            if (coeff > 0.)
                                            {
                                                if (isActive(selfLabel))
                                                {
                                                    localActiveDivergenceMatrixElements.push_back(Eigen::Triplet<SolveReal>(selfActiveIndex, edgeStressIndex, contribution));

                                                    // solid boundary
                                                    if (getFieldValue(*edgeFluidWeights(edgeAxis), edge) < 1.)
                                                    {
                                                        SolveReal solidEdgeVolume = 1. - getFieldValue(*edgeFluidWeights(edgeAxis), edge);
                                                        SolveReal liquidEdgeVolume = getFieldValue(*edgeLiquidWeights(edgeAxis), edge);
                                                        SolveReal solidCoeff = liquidEdgeVolume * solidEdgeVolume * invDx;
                                                        SolveReal solidContribution = divergenceSign * coeff;
                                                        SolveReal svel = getFieldValue(*myCollisionVelocityField->getField(faceAxis), face);
                                                        localStressRHSElements.push_back(Eigen::Triplet<SolveReal>(edgeStressIndex, 0, -1. * solidContribution * svel));
                                                    }
                                                    if (getFieldValue(*faceFluidWeights(faceAxis), face) < 1.)
                                                    { // RHS Contribution
                                                        SolveReal solidFaceVolume = 1. - getFieldValue(*faceFluidWeights(faceAxis), face);
                                                        SolveReal liquidEdgeVolume = getFieldValue(*edgeLiquidWeights(edgeAxis), edge);
                                                        SolveReal solidCoeff = liquidEdgeVolume * solidFaceVolume * invDx;
                                                        SolveReal solidContribution = divergenceSign * coeff;
                                                        SolveReal svel = getFieldValue(*myCollisionVelocityField->getField(faceAxis), face);
                                                        localStressRHSElements.push_back(Eigen::Triplet<SolveReal>(edgeStressIndex, 0, solidContribution * svel));
                                                    }
                                                }
                                                else if (isReduced(selfLabel))
                                                {
                                                    // setup offsets
                                                    UT_Vector3T<SolveReal> selfOffset;
                                                    UT_Vector3T<SolveReal> facePos(face);
                                                    facePos[faceAxis] -= 0.5;
                                                    facePos *= dx;
                                                    selfOffset = facePos - reducedRegionCOM[selfReducedIndex];

                                                    ColumnVector colVec = buildConversionCoefficients(selfOffset, faceAxis);

                                                    for (int n = 0; n < REDUCED_DOF; n++)
                                                        localReducedDivergenceMatrixElements.push_back(Eigen::Triplet<SolveReal>(REDUCED_DOF * selfReducedIndex + n, edgeStressIndex, contribution * colVec(n)));
                                                }
                                            }
                                        }
                                        if (isReduced(edgeLabel))
                                        {
                                            exint reducedEdgeIndex = getFieldValue(*edgeReducedIndices(edgeAxis), edge);
                                            exint reducedEdgeStressIndex = reducedEdgeStressDOF(reducedEdgeIndex, edgeAxis);
                                            SolveReal contribution = -1. * invDx * divergenceSign;

                                            if (isReduced(selfLabel))
                                            {
                                                // setup offsets
                                                UT_Vector3T<SolveReal> selfOffset;
                                                UT_Vector3T<SolveReal> facePos(face);
                                                facePos[faceAxis] -= 0.5;
                                                facePos *= dx;
                                                selfOffset = facePos - reducedRegionCOM[selfReducedIndex];

                                                ColumnVector colVec = buildConversionCoefficients(selfOffset, faceAxis);

                                                for (int n = 0; n < REDUCED_DOF; n++)
                                                    localReducedInternalDivergenceMatrixElements.push_back(Eigen::Triplet<SolveReal>(REDUCED_DOF * selfReducedIndex + n, reducedEdgeStressIndex, contribution * colVec(n)));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                return 0;
            }
        );
    }

    // setup edge stress diagonal terms
    for (int edgeAxis : {0, 1, 2}) {

        UT_Interrupt* boss = UTgetInterrupt();

        UT_ThreadedAlgorithm buildLinearSystemAlgorithm;
        buildLinearSystemAlgorithm.run([&](const UT_JobInfo& info)
            {
                Triplets& localActiveStressMatrixElements = parallelActiveStressMatrixElements[info.job()];
                Triplets& localActiveStressInverseMatrixElements = parallelActiveStressInverseMatrixElements[info.job()];
                Triplets& localReducedStressMatrixElements = parallelReducedStressMatrixElements[info.job()];
                Triplets& localReducedStressInverseMatrixElements = parallelReducedStressInverseMatrixElements[info.job()];

                UT_VoxelArrayIteratorI vit;
                vit.setConstArray(edgeLabels(edgeAxis)->field());
                vit.splitByTile(info);

                UT_VoxelTileIteratorI vitt;
                for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
                {
                    if (boss->opInterrupt())
                        break;

                    if (!vit.isTileConstant() || isSolved(vit.getValue()))
                    {
                        vitt.setTile(vit);
                        for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                        {
                            int j = edgeAxis;
                            UT_Vector3I edge(vitt.x(), vitt.y(), vitt.z());

                            // we need to make sure we translate the axis DOF into the index of the entire stress DOFs
                            MaterialLabels edgeLabel = (MaterialLabels)getFieldValue(*edgeLabels(edgeAxis), edge);
                            exint edgeStressIndex = edgeStressDOF(getFieldValue(*edgeActiveIndices(edgeAxis), edge), edgeAxis);

                            // grab volume weight
                            // in progress clamping for fluid weight
                            SolveReal volumeWeight =
                                //(SolveReal)getFieldValue(*edgeFluidWeights(edgeAxis), edge)
                                SYSclamp((SolveReal)getFieldValue(*edgeFluidWeights(edgeAxis), edge), MINWEIGHT, 1.0)
                                * (SolveReal)getFieldValue(*edgeLiquidWeights(edgeAxis), edge);

                            // grab viscosity
                            UT_Vector3 point;
                            edgeLabels(edgeAxis)->indexToPos(edge[0], edge[1], edge[2], point);
                            SolveReal localViscosity = getLocalViscosity(point);

                            SolveReal invLocalViscosity = SYSclamp(1. / localViscosity, 0., 1e10);

                            // setup diagonal stress identity matrix
                            if (isActive(edgeLabel))
                            {
                                localActiveStressInverseMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    edgeStressIndex,
                                    edgeStressIndex,
                                    2. * invLocalViscosity * volumeWeight)
                                );
                                localActiveStressMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    edgeStressIndex,
                                    edgeStressIndex,
                                    0.5 * localViscosity * SYSclamp(1. / volumeWeight, 0., 1.e2))
                                );
                            }
                            if (isReduced(edgeLabel))
                            {
                                exint reducedEdgeIndex = getFieldValue(*edgeReducedIndices(edgeAxis), edge);
                                exint reducedEdgeStressIndex = reducedEdgeStressDOF(reducedEdgeIndex, edgeAxis);

                                localReducedStressInverseMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    reducedEdgeStressIndex,
                                    reducedEdgeStressIndex,
                                    2. * invLocalViscosity * volumeWeight)
                                );
                                localReducedStressMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    reducedEdgeStressIndex,
                                    reducedEdgeStressIndex,
                                    0.5 * localViscosity * SYSclamp(1. / volumeWeight, 0., 1.e2))
                                );
                            }
                        }
                    }
                }

                return 0;
            });
    }
    // setup center stress diagonal terms
    {
        UT_Interrupt* boss = UTgetInterrupt();

        UT_ThreadedAlgorithm buildLinearSystemAlgorithm;
        buildLinearSystemAlgorithm.run([&](const UT_JobInfo& info)
            {
                Triplets& localActiveStressMatrixElements = parallelActiveStressMatrixElements[info.job()];
                Triplets& localActiveStressInverseMatrixElements = parallelActiveStressInverseMatrixElements[info.job()];
                Triplets& localReducedStressMatrixElements = parallelReducedStressMatrixElements[info.job()];
                Triplets& localReducedStressInverseMatrixElements = parallelReducedStressInverseMatrixElements[info.job()];

                UT_VoxelArrayIteratorI vit;
                vit.setConstArray(centerLabels.field());
                vit.splitByTile(info);

                UT_VoxelTileIteratorI vitt;
                for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
                {
                    if (boss->opInterrupt())
                        break;

                    if (!vit.isTileConstant() || isSolved(vit.getValue()))
                    {
                        vitt.setTile(vit);
                        for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                        {
                            UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

                            // we need to make sure we translate the axis DOF into the index of the entire stress DOFs
                            MaterialLabels cellLabel = (MaterialLabels)getFieldValue(centerLabels, cell);
                            exint cellIncompleteIndex = getFieldValue(centerActiveIndices, cell);
                            exint cellStressXXIndex = centerStressDOF(cellIncompleteIndex, StressType::XX);
                            exint cellStressYYIndex = centerStressDOF(cellIncompleteIndex, StressType::YY);
                            exint cellStressZZIndex = centerStressDOF(cellIncompleteIndex, StressType::ZZ);

                            // grab volume weight
                            // in progress clamping
                            SolveReal volumeWeight =
                                //(SolveReal)getFieldValue(centerFluidWeights, cell)
                                SYSclamp((SolveReal)getFieldValue(centerFluidWeights, cell), MINWEIGHT, 1.0)
                                * (SolveReal)getFieldValue(centerLiquidWeights, cell);

                            // grab viscosity
                            UT_Vector3 point;
                            centerLabels.indexToPos(cell[0], cell[1], cell[2], point);
                            SolveReal localViscosity = getLocalViscosity(point);

                            SolveReal invLocalViscosity = SYSclamp(1. / localViscosity, 0., 1.e10);

                            // setup diagonal mass matrix, rhs, and warm start guess
                            if (isActive(cellLabel))
                            {
                                localActiveStressInverseMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    cellStressXXIndex,
                                    cellStressXXIndex,
                                    invLocalViscosity * SYSclamp(volumeWeight, 1.e-2, 1.))
                                );
                                localActiveStressInverseMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    cellStressYYIndex,
                                    cellStressYYIndex,
                                    invLocalViscosity * SYSclamp(volumeWeight, 1.e-2, 1.))
                                );
                                localActiveStressInverseMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    cellStressZZIndex,
                                    cellStressZZIndex,
                                    invLocalViscosity * SYSclamp(volumeWeight, 1.e-2, 1.))
                                );

                                localActiveStressMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    cellStressXXIndex,
                                    cellStressXXIndex,
                                    localViscosity * SYSclamp(1. / volumeWeight, 0., 1.e2))
                                );
                                localActiveStressMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    cellStressYYIndex,
                                    cellStressYYIndex,
                                    localViscosity * SYSclamp(1. / volumeWeight, 0., 1.e2))
                                );
                                localActiveStressMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    cellStressZZIndex,
                                    cellStressZZIndex,
                                    localViscosity * SYSclamp(1. / volumeWeight, 0., 1.e2))
                                );
                            }
                            if (isReduced(cellLabel))
                            {
                                exint reducedCellIndex = getFieldValue(centerReducedIndices, cell);
                                exint reducedCellStressXXIndex = reducedCenterStressDOF(reducedCellIndex, StressType::XX);
                                exint reducedCellStressYYIndex = reducedCenterStressDOF(reducedCellIndex, StressType::YY);
                                exint reducedCellStressZZIndex = reducedCenterStressDOF(reducedCellIndex, StressType::ZZ);

                                localReducedStressInverseMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    reducedCellStressXXIndex,
                                    reducedCellStressXXIndex,
                                    invLocalViscosity * SYSclamp(volumeWeight, 1.e-2, 1.))
                                );
                                localReducedStressInverseMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    reducedCellStressYYIndex,
                                    reducedCellStressYYIndex,
                                    invLocalViscosity * SYSclamp(volumeWeight, 1.e-2, 1.))
                                );
                                localReducedStressInverseMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    reducedCellStressZZIndex,
                                    reducedCellStressZZIndex,
                                    invLocalViscosity * SYSclamp(volumeWeight, 1.e-2, 1.))
                                );

                                localReducedStressMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    reducedCellStressXXIndex,
                                    reducedCellStressXXIndex,
                                    localViscosity * SYSclamp(1. / volumeWeight, 0., 1.e2))
                                );
                                localReducedStressMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    reducedCellStressYYIndex,
                                    reducedCellStressYYIndex,
                                    localViscosity * SYSclamp(1. / volumeWeight, 0., 1.e2))
                                );
                                localReducedStressMatrixElements.push_back(Eigen::Triplet<SolveReal>(
                                    reducedCellStressZZIndex,
                                    reducedCellStressZZIndex,
                                    localViscosity * SYSclamp(1. / volumeWeight, 0., 1.e2))
                                );
                            }
                        }
                    }
                }

                return 0;
            });

    }
}