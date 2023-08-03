#include "HDK_PolyStokesSolver.h"
//#include "CuConjugateGradient.h"

// for constructing reduced indices
//#include <SIM/SIM_VolumetricConnectedComponentBuilder.h>
// for printing out data
#include <SIM/SIM_GeometryCopy.h>
// saving out matrix
#include <unsupported/Eigen/SparseExtra>
// solvers
#include "pcg.h"
// cumat
//#include <cuMat/Core>
//#include <cuMat/IterativeLinearSolvers>
//#include <cuMat/Sparse>

// Constructor and destructor
HDK_PolyStokes::Solver::Solver(
    HDK_PolyStokes& _parent,
    const fpreal _dx, const fpreal _dt,
    SIM_VectorField* _velocityField, 
    const SIM_VectorField* _collisionVelocityField,
    const SIM_ScalarField* _surfaceField,
    const SIM_VectorField* _surfaceWeights,
    const SIM_ScalarField* _collisionField,
    const SIM_VectorField* _collisionWeights,
    const SIM_ScalarField* _densityField,
    const fpreal _constantDensity,
    const SIM_ScalarField* _viscosityField
) :
    myParent(_parent)
    , myVelocityField(_velocityField)
    , myCollisionVelocityField(_collisionVelocityField)
    , dx(_dx)
    , invDx(1./dx)
    , invDx2(1./(dx*dx))
    , dt(_dt)
    , invDt(1./dt)
    , size(_velocityField->getSize())
    , orig(_velocityField->getOrig())
    , resolution(_velocityField->getTotalVoxelRes())
    , nx(resolution.x())
    , ny(resolution.y())
    , nz(resolution.z())
    , myUseSurfaceWeights(_parent.getUseInputSurfaceWeights())
    , mySurfaceField(_surfaceField)
    , mySurfaceWeights(_surfaceWeights)
    , myUseCollisionWeights(_parent.getUseInputCollisionWeights())
    , myCollisionField(_collisionField)
    , myCollisionWeights(_collisionWeights)
    , mySurfaceFieldData(_surfaceField->getField())
    , myCollisionFieldData(_collisionField->getField())
    , myViscosityFieldData(_viscosityField->getField())
    , myDensityFieldData(_densityField->getField())
    , myMinDensity(_parent.getMinDensity())
    , myMaxDensity(_parent.getMaxDensity())
    , myConstantDensity(_constantDensity)
    , myLiquidBoundaryLayerSize(_parent.getActiveLiquidBoundaryLayerSize())
    , mySolidBoundaryLayerSize(_parent.getActiveSolidBoundaryLayerSize())
    , mySolverTolerance(_parent.getSolverTolerance())
    , mySolverMaxIterations(_parent.getSolverMaxIterations())
    , myDoReducedRegions(_parent.getDoReducedRegions())
    , myDoTile(_parent.getDoTile())
    , myTileSize(_parent.getTileSize())
    , myTilePadding(_parent.getTilePadding())
    , myMatrixScheme(_parent.getMatrixScheme())
    , mySolverType(_parent.getSolverType())
{
    // create and initialize fields
    // weights
    initField(centerLiquidWeights, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::CENTER);
    initField(edgeXYLiquidWeights, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::EDGEXY);
    initField(edgeXZLiquidWeights, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::EDGEXZ);
    initField(edgeYZLiquidWeights, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::EDGEYZ);
    initField(faceXLiquidWeights, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::FACEX);
    initField(faceYLiquidWeights, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::FACEY);
    initField(faceZLiquidWeights, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::FACEZ);
    initField(centerFluidWeights, *myCollisionFieldData, (SIM_FieldSample)SamplingType::CENTER);
    initField(edgeXYFluidWeights, *myCollisionFieldData, (SIM_FieldSample)SamplingType::EDGEXY);
    initField(edgeXZFluidWeights, *myCollisionFieldData, (SIM_FieldSample)SamplingType::EDGEXZ);
    initField(edgeYZFluidWeights, *myCollisionFieldData, (SIM_FieldSample)SamplingType::EDGEYZ);
    initField(faceXFluidWeights, *myCollisionFieldData, (SIM_FieldSample)SamplingType::FACEX);
    initField(faceYFluidWeights, *myCollisionFieldData, (SIM_FieldSample)SamplingType::FACEY);
    initField(faceZFluidWeights, *myCollisionFieldData, (SIM_FieldSample)SamplingType::FACEZ);

    // labels
    initField(centerLabels, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::CENTER);
    initField(edgeXYLabels, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::EDGEXY);
    initField(edgeXZLabels, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::EDGEXZ);
    initField(edgeYZLabels, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::EDGEYZ);
    initField(faceXLabels, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::FACEX);
    initField(faceYLabels, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::FACEY);
    initField(faceZLabels, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::FACEZ);
    centerLabels.makeConstant(MaterialLabels::UNASSIGNED);
    edgeXYLabels.makeConstant(MaterialLabels::UNASSIGNED);
    edgeXZLabels.makeConstant(MaterialLabels::UNASSIGNED);
    edgeYZLabels.makeConstant(MaterialLabels::UNASSIGNED);
    faceXLabels.makeConstant(MaterialLabels::UNASSIGNED);
    faceYLabels.makeConstant(MaterialLabels::UNASSIGNED);
    faceZLabels.makeConstant(MaterialLabels::UNASSIGNED);
    centerLabels.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    edgeXYLabels.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    edgeXZLabels.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    edgeYZLabels.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    faceXLabels.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    faceYLabels.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    faceZLabels.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);

    //centerReducedIndices.match(*mySurfaceFieldData);
    initField(centerReducedIndices, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::CENTER);
    initField(edgeXYReducedIndices, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::EDGEXY);
    initField(edgeXZReducedIndices, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::EDGEXZ);
    initField(edgeYZReducedIndices, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::EDGEYZ);
    initField(faceXReducedIndices, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::FACEX);
    initField(faceYReducedIndices, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::FACEY);
    initField(faceZReducedIndices, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::FACEZ);
    centerReducedIndices.makeConstant(MaterialLabels::UNASSIGNED);
    edgeXYReducedIndices.makeConstant(MaterialLabels::UNASSIGNED);
    edgeXZReducedIndices.makeConstant(MaterialLabels::UNASSIGNED);
    edgeYZReducedIndices.makeConstant(MaterialLabels::UNASSIGNED);
    faceXReducedIndices.makeConstant(MaterialLabels::UNASSIGNED);
    faceYReducedIndices.makeConstant(MaterialLabels::UNASSIGNED);
    faceZReducedIndices.makeConstant(MaterialLabels::UNASSIGNED);
    centerReducedIndices.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    edgeXYReducedIndices.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    edgeXZReducedIndices.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    edgeYZReducedIndices.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    faceXReducedIndices.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    faceYReducedIndices.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    faceZReducedIndices.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);

    initField(centerActiveIndices, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::CENTER);
    initField(edgeXYActiveIndices, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::EDGEXY);
    initField(edgeXZActiveIndices, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::EDGEXZ);
    initField(edgeYZActiveIndices, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::EDGEYZ);
    initField(faceXActiveIndices, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::FACEX);
    initField(faceYActiveIndices, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::FACEY);
    initField(faceZActiveIndices, *mySurfaceFieldData, (SIM_FieldSample)SamplingType::FACEZ);
    centerActiveIndices.makeConstant(MaterialLabels::UNASSIGNED);
    edgeXYActiveIndices.makeConstant(MaterialLabels::UNASSIGNED);
    edgeXZActiveIndices.makeConstant(MaterialLabels::UNASSIGNED);
    edgeYZActiveIndices.makeConstant(MaterialLabels::UNASSIGNED);
    faceXActiveIndices.makeConstant(MaterialLabels::UNASSIGNED);
    faceYActiveIndices.makeConstant(MaterialLabels::UNASSIGNED);
    faceZActiveIndices.makeConstant(MaterialLabels::UNASSIGNED);
    centerActiveIndices.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    edgeXYActiveIndices.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    edgeXZActiveIndices.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    edgeYZActiveIndices.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    faceXActiveIndices.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    faceYActiveIndices.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    faceZActiveIndices.setBorder(UT_VOXELBORDER_CONSTANT, MaterialLabels::UNASSIGNED);
    // create and initialize multithreading vars
    myThreadCount = UT_Thread::getNumProcessors();
}

HDK_PolyStokes::Solver::~Solver()
{
}

void
HDK_PolyStokes::Solver::buildIntegrationWeights()
{
    // liquid weights -> 1 for inside liquid and 0 for air
    liquidWeights = {
        &centerLiquidWeights,
        &edgeXYLiquidWeights,
        &edgeXZLiquidWeights,
        &edgeYZLiquidWeights,
        NULL, NULL, NULL
    };
    // fluid weights -> 1 for inside fluid and 0 for solid
    fluidWeights = {
        &centerFluidWeights,
        &edgeXYFluidWeights,
        &edgeXZFluidWeights,
        &edgeYZFluidWeights,
        NULL, NULL, NULL
    };

    // compute liquid/air weights
    // liquid = 1, outside = 0
    if (myUseSurfaceWeights)
    {
        // put in the input face weights
        faceXLiquidWeights = *mySurfaceWeights->getField(0);
        faceYLiquidWeights = *mySurfaceWeights->getField(1);
        faceZLiquidWeights = *mySurfaceWeights->getField(2);
        liquidWeights[4] = &faceXLiquidWeights;
        liquidWeights[5] = &faceYLiquidWeights;
        liquidWeights[6] = &faceZLiquidWeights;
        for (int i = 4; i < 7; i++)
            liquidWeights[i]->setScaleDivideThreshold(1., NULL, NULL, MINWEIGHT);
    }
    else
    {
        faceXLiquidWeights.computeSDFWeightsSampled(mySurfaceFieldData, NSAMPLES, false, MINWEIGHT);
        faceYLiquidWeights.computeSDFWeightsSampled(mySurfaceFieldData, NSAMPLES, false, MINWEIGHT);
        faceZLiquidWeights.computeSDFWeightsSampled(mySurfaceFieldData, NSAMPLES, false, MINWEIGHT);
        liquidWeights[4] = &faceXLiquidWeights;
        liquidWeights[5] = &faceYLiquidWeights;
        liquidWeights[6] = &faceZLiquidWeights;
    }
    centerLiquidWeights.computeSDFWeightsSampled(mySurfaceFieldData, NSAMPLES, false, MINWEIGHT);
    edgeXYLiquidWeights.computeSDFWeightsSampled(mySurfaceFieldData, NSAMPLES, false, MINWEIGHT);
    edgeXZLiquidWeights.computeSDFWeightsSampled(mySurfaceFieldData, NSAMPLES, false, MINWEIGHT);
    edgeYZLiquidWeights.computeSDFWeightsSampled(mySurfaceFieldData, NSAMPLES, false, MINWEIGHT);

    // compute fluid/solid weights
    // fluid = 1, solid = 0
    if (myUseCollisionWeights)
    {
        // put in the input face weights
        faceXFluidWeights = *myCollisionWeights->getField(0);
        faceYFluidWeights = *myCollisionWeights->getField(1);
        faceZFluidWeights = *myCollisionWeights->getField(2);
        fluidWeights[4] = &faceXFluidWeights;
        fluidWeights[5] = &faceYFluidWeights;
        fluidWeights[6] = &faceZFluidWeights;
        for (int i = 4; i < 7; i++)
            fluidWeights[i]->setScaleDivideThreshold(1., NULL, NULL, MINWEIGHT);
    }
    else
    {
        faceXFluidWeights.computeSDFWeightsSampled(myCollisionFieldData, NSAMPLES, false, 0.);
        faceYFluidWeights.computeSDFWeightsSampled(myCollisionFieldData, NSAMPLES, false, 0.);
        faceZFluidWeights.computeSDFWeightsSampled(myCollisionFieldData, NSAMPLES, false, 0.);
        fluidWeights[4] = &faceXFluidWeights;
        fluidWeights[5] = &faceYFluidWeights;
        fluidWeights[6] = &faceZFluidWeights;
    }
    centerFluidWeights.computeSDFWeightsSampled(myCollisionFieldData, NSAMPLES, false, 0.);
    edgeXYFluidWeights.computeSDFWeightsSampled(myCollisionFieldData, NSAMPLES, false, 0.);
    edgeXZFluidWeights.computeSDFWeightsSampled(myCollisionFieldData, NSAMPLES, false, 0.);
    edgeYZFluidWeights.computeSDFWeightsSampled(myCollisionFieldData, NSAMPLES, false, 0.);
}

void
HDK_PolyStokes::Solver::buildIntegrationWeightsAlt()
{
    // todo make this a tunable param
    const int numberOfSamples = 2;
    const SolveReal extrapolation = dx * 0.;

    {
        computeIntegrationWeights(centerLiquidWeights, *mySurfaceFieldData, SIM_SAMPLE_CENTER, numberOfSamples);

        computeIntegrationWeights(faceXLiquidWeights, *mySurfaceFieldData, SIM_SAMPLE_FACEX, numberOfSamples);
        computeIntegrationWeights(faceYLiquidWeights, *mySurfaceFieldData, SIM_SAMPLE_FACEY, numberOfSamples);
        computeIntegrationWeights(faceZLiquidWeights, *mySurfaceFieldData, SIM_SAMPLE_FACEZ, numberOfSamples);

        computeIntegrationWeights(edgeYZLiquidWeights, *mySurfaceFieldData, SIM_SAMPLE_EDGEYZ, numberOfSamples);
        computeIntegrationWeights(edgeXZLiquidWeights, *mySurfaceFieldData, SIM_SAMPLE_EDGEXZ, numberOfSamples);
        computeIntegrationWeights(edgeXYLiquidWeights, *mySurfaceFieldData, SIM_SAMPLE_EDGEXY, numberOfSamples);
    }

    {
        computeSolidIntegrationWeights(centerFluidWeights, *mySurfaceFieldData, *myCollisionFieldData, SIM_SAMPLE_CENTER, numberOfSamples, -extrapolation);

        computeSolidIntegrationWeights(faceXFluidWeights, *mySurfaceFieldData, *myCollisionFieldData, SIM_SAMPLE_FACEX, numberOfSamples, -extrapolation);
        computeSolidIntegrationWeights(faceYFluidWeights, *mySurfaceFieldData, *myCollisionFieldData, SIM_SAMPLE_FACEY, numberOfSamples, -extrapolation);
        computeSolidIntegrationWeights(faceZFluidWeights, *mySurfaceFieldData, *myCollisionFieldData, SIM_SAMPLE_FACEZ, numberOfSamples, -extrapolation);

        computeSolidIntegrationWeights(edgeYZFluidWeights, *mySurfaceFieldData, *myCollisionFieldData, SIM_SAMPLE_EDGEYZ, numberOfSamples, -extrapolation);
        computeSolidIntegrationWeights(edgeXZFluidWeights, *mySurfaceFieldData, *myCollisionFieldData, SIM_SAMPLE_EDGEXZ, numberOfSamples, -extrapolation);
        computeSolidIntegrationWeights(edgeXYFluidWeights, *mySurfaceFieldData, *myCollisionFieldData, SIM_SAMPLE_EDGEXY, numberOfSamples, -extrapolation);
    }

    // liquid weights -> 1 for inside liquid and 0 for air
    liquidWeights = {
        &centerLiquidWeights,
        &edgeXYLiquidWeights,
        &edgeXZLiquidWeights,
        &edgeYZLiquidWeights,
        &faceXLiquidWeights, 
        &faceYLiquidWeights,
        &faceZLiquidWeights
    };
    // fluid weights -> 1 for inside fluid and 0 for solid
    fluidWeights = {
        &centerFluidWeights,
        &edgeXYFluidWeights,
        &edgeXZFluidWeights,
        &edgeYZFluidWeights,
        &faceXFluidWeights,
        &faceYFluidWeights,
        &faceZFluidWeights
    };
}

void
HDK_PolyStokes::Solver::computeIntegrationWeights(SIM_RawField& integrationWeights,
    const SIM_RawField& liquidSurface,
    const SIM_FieldSample sample,
    const int numberOfSamples)
{
    UT_Vector3 size = liquidSurface.getSize(), orig = liquidSurface.getOrig();

    UT_Vector3i voxelRes;
    liquidSurface.getVoxelRes(voxelRes[0], voxelRes[1], voxelRes[2]);

    integrationWeights.init(sample, orig, size, voxelRes[0], voxelRes[1], voxelRes[2]);
    integrationWeights.makeConstant(0);
    integrationWeights.computeSDFWeightsSampled(&liquidSurface, numberOfSamples, false, 0);
}

void
HDK_PolyStokes::Solver::computeSolidIntegrationWeights(SIM_RawField& solidIntegrationWeights,
    const SIM_RawField& liquidSurface,
    const SIM_RawField& solidSurface,
    const SIM_FieldSample sample,
    const int numberOfSamples,
    const fpreal extrapolation)
{
    UT_Vector3 size = liquidSurface.getSize(), orig = liquidSurface.getOrig();

    UT_Vector3i voxelRes;
    liquidSurface.getVoxelRes(voxelRes[0], voxelRes[1], voxelRes[2]);

    solidIntegrationWeights.init(sample, orig, size, voxelRes[0], voxelRes[1], voxelRes[2]);
    solidIntegrationWeights.makeConstant(0);
    solidIntegrationWeights.computeSDFWeightsSampled(&solidSurface, numberOfSamples,
        false, 0, -extrapolation);
    // TOCHECK: last param is negative or positive
    // TOCHECK: 3rd param is invert in case we need to invert the solid weights
}

void
HDK_PolyStokes::Solver::computeCenterOfMasses()
{
    reducedRegionCOM.setSize(myInteriorRegionCount);
    reducedRegionCOM.constant(UT_Vector3T<SolveReal>(0, 0, 0));

    UT_Array<UT_Array<UT_Vector3T<SolveReal>>> parallelInteriorRegionCOM;
    parallelInteriorRegionCOM.setSize(myThreadCount);

    UT_Array<UT_Array<SolveReal>> parallelInteriorRegionCellCount;
    parallelInteriorRegionCellCount.setSize(myThreadCount);

    for (int thread = 0; thread < myThreadCount; ++thread)
    {
        parallelInteriorRegionCOM[thread].setSize(myInteriorRegionCount);
        parallelInteriorRegionCOM[thread].constant(UT_Vector3T<SolveReal>(0, 0, 0));

        parallelInteriorRegionCellCount[thread].setSize(myInteriorRegionCount);
        parallelInteriorRegionCellCount[thread].constant(0);
    }

    buildReducedRegionCOM(parallelInteriorRegionCOM, parallelInteriorRegionCellCount);

    UT_Array<SolveReal> interiorRegionCellCount;
    interiorRegionCellCount.setSize(myInteriorRegionCount);
    interiorRegionCellCount.constant(0);

    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
            {
                for (int thread = 0; thread < myThreadCount; ++thread)
                {
                    reducedRegionCOM[interiorRegion] += parallelInteriorRegionCOM[thread][interiorRegion];
                    interiorRegionCellCount[interiorRegion] += parallelInteriorRegionCellCount[thread][interiorRegion];
                }
            }
        });

    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
                reducedRegionCOM[interiorRegion] *= dx / interiorRegionCellCount[interiorRegion];
        });
}

void
HDK_PolyStokes::Solver::computeLeastSquaresFits()
{
    UT_Array<UT_Array<ReducedMatrix>> parallelBestFitMatrix;
    parallelBestFitMatrix.setSize(myThreadCount);
    UT_Array<UT_Array<ColumnVector>> parallelBestFitRHS;
    parallelBestFitRHS.setSize(myThreadCount);
    for (int thread = 0; thread < myThreadCount; ++thread)
    {
        parallelBestFitMatrix[thread].setSize(myInteriorRegionCount);
        parallelBestFitMatrix[thread].constant(ReducedMatrix::Zero());
        parallelBestFitRHS[thread].setSize(myInteriorRegionCount);
        parallelBestFitRHS[thread].constant(ColumnVector::Zero());
    }

    buildInteriorBestFitSystems(parallelBestFitMatrix, parallelBestFitRHS);

    reducedRegionBestFitVectors.setSize(myInteriorRegionCount);
    reducedRegionBestFitVectors.constant(ColumnVector::Zero());
    // Compile best fit matrices and rhs vectors
    UT_Array<ReducedMatrix> interiorBestFitMatrix;
    interiorBestFitMatrix.setSize(myInteriorRegionCount);
    interiorBestFitMatrix.constant(ReducedMatrix::Zero());
    UT_Array<ColumnVector> interiorBestFitRHS;
    interiorBestFitRHS.setSize(myInteriorRegionCount);
    interiorBestFitRHS.constant(ColumnVector::Zero());

    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
            {
                for (int thread = 0; thread < myThreadCount; ++thread)
                {
                    interiorBestFitMatrix[interiorRegion] += parallelBestFitMatrix[thread][interiorRegion];
                    interiorBestFitRHS[interiorRegion] += parallelBestFitRHS[thread][interiorRegion];
                }
            }
        });
    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
                reducedRegionBestFitVectors[interiorRegion] = interiorBestFitMatrix[interiorRegion].fullPivLu().solve(interiorBestFitRHS[interiorRegion]);
        });
}

void
HDK_PolyStokes::Solver::computeReducedMassMatrices()
{
    UT_Array<UT_Array<ReducedMatrix>> parallelReducedMassMatrix;
    parallelReducedMassMatrix.setSize(myThreadCount);
    for (int thread = 0; thread < myThreadCount; ++thread)
    {
        parallelReducedMassMatrix[thread].setSize(myInteriorRegionCount);
        parallelReducedMassMatrix[thread].constant(ReducedMatrix::Zero());
    }

    buildReducedMassMatrixSystems(parallelReducedMassMatrix);
    // Compile interior mass matrix elements
    reducedMassMatrices.setSize(myInteriorRegionCount);
    reducedMassMatrices.constant(ReducedMatrix::Zero());

    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
                for (int thread = 0; thread < myThreadCount; ++thread)
                    reducedMassMatrices[interiorRegion] += parallelReducedMassMatrix[thread][interiorRegion];
        });
}

void
HDK_PolyStokes::Solver::computeReducedViscosityMatrices()
{
    UT_Array<UT_Array<ReducedMatrix>> parallelReducedViscosityMatrix;
    parallelReducedViscosityMatrix.setSize(myThreadCount);

    for (int thread = 0; thread < myThreadCount; ++thread)
    {
        parallelReducedViscosityMatrix[thread].setSize(myInteriorRegionCount);
        parallelReducedViscosityMatrix[thread].constant(ReducedMatrix::Zero());
    }

    buildReducedViscosityMatrixSystems(parallelReducedViscosityMatrix);
    // Compile interior viscosity matrix elements
    reducedViscosityMatrices.setSize(myInteriorRegionCount);
    reducedViscosityMatrices.constant(ReducedMatrix::Zero());

    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
                for (int thread = 0; thread < myThreadCount; ++thread)
                    reducedViscosityMatrices[interiorRegion] += parallelReducedViscosityMatrix[thread][interiorRegion];
        });
}

void
HDK_PolyStokes::Solver::computeReducedViscosityMatricesInteriorOnly()
{
    UT_Array<UT_Array<ReducedMatrix>> parallelReducedViscosityMatrix;
    parallelReducedViscosityMatrix.setSize(myThreadCount);
    for (int thread = 0; thread < myThreadCount; ++thread)
    {
        parallelReducedViscosityMatrix[thread].setSize(myInteriorRegionCount);
        parallelReducedViscosityMatrix[thread].constant(ReducedMatrix::Zero());
    }

    buildReducedViscosityMatrixSystemsInteriorOnly(parallelReducedViscosityMatrix);
    // Compile interior viscosity matrix elements
    reducedViscosityMatrices.setSize(myInteriorRegionCount);
    reducedViscosityMatrices.constant(ReducedMatrix::Zero());

    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
                for (int thread = 0; thread < myThreadCount; ++thread)
                    reducedViscosityMatrices[interiorRegion] += parallelReducedViscosityMatrix[thread][interiorRegion];
        });
}

void
HDK_PolyStokes::Solver::recoverVelocityFromPressureStress()
{
    exint pressureOffset = 0;
    exint stressOffset = nPressures;
    // todo introduce a check here to make sure we're in the correct mode (ie nSystemSize matches what is expected)
    nSystemSize = nPressures + nStresses;
    exint velocitySystemSize = nActiveVs + nReducedVs;

    Vector pressureStressSolution = solutionVector;
    Vector pressureSolution = pressureStressSolution.head(nPressures);
    Vector stressSolution = pressureStressSolution.tail(nStresses);

    // right now, applySolutionToVelocity expects the velocities to be at the top of the solution vector, so use this
    solutionVector = Vector::Zero(velocitySystemSize);
    solutionVector.head(nActiveVs) = dt * McInv_Matrix * (invDt * activeRHSVector - G_Matrix * pressureSolution - Dt_Matrix * stressSolution);
    // todo
    solutionVector.tail(nReducedVs) = Inv_Mr_plus_2JDtuDJ_Matrix * (invDt * reducedRHSVector - JG_Matrix * pressureSolution - JDt_Matrix * stressSolution);
}

void
HDK_PolyStokes::Solver::initializeGuessVectors()
{
    activeGuessVector   = Vector::Zero(nActiveVs);
    reducedGuessVector  = Vector::Zero(nReducedVs);
    pressureGuessVector = Vector::Zero(nPressures);
    stressGuessVector   = Vector::Zero(nStresses);
}

void
HDK_PolyStokes::Solver::constructGuessVectors()
{
    for (exint i = 0; i < nActiveVs; ++i)
        activeGuessVector(i) = oldActiveVs.coeffRef(i, 0);
    for (exint i = 0; i < myInteriorRegionCount; ++i)
        for (exint j = 0; j < REDUCED_DOF; ++j)
            reducedGuessVector(i * REDUCED_DOF + j) = reducedRegionBestFitVectors[i](j);
    pressureGuessVector = -G_Matrix.transpose() * activeGuessVector - JG_Matrix.transpose() * reducedGuessVector;
    stressGuessVector = -2. * uInv_Matrix * (-Dt_Matrix.transpose() * activeGuessVector - JDt_Matrix.transpose() * reducedGuessVector);
}

void
HDK_PolyStokes::Solver::exportMatrices(std::string prefix)
{
    // write out matrices for debugging
    Eigen::saveMarket(A, prefix + "Mat_A.mtx");
    //Eigen::saveMarket(preconditioner, prefix + "Mat_preconditioner.mtx"); // we got rid of explicit matrix preconditioners, need to call the class now
    Eigen::saveMarketVector(b, prefix + "Vec_b.mtx");
    Eigen::saveMarketVector(guessVector, prefix + "Vec_guess.mtx");
}

void
HDK_PolyStokes::Solver::exportComponentMatrices(std::string prefix)
{
    Eigen::saveMarket(Mc_Matrix, prefix + "Mat_Mc.mtx");
    Eigen::saveMarket(McInv_Matrix, prefix + "Mat_McInv.mtx");
    Eigen::saveMarket(Mr_Matrix, prefix + "Mat_Mr.mtx");
    Eigen::saveMarket(MrInv_Matrix, prefix + "Mat_MrInv.mtx");
    Eigen::saveMarket(JDtuDJ_Matrix, prefix + "Mat_JDtuDJ.mtx");
    Eigen::saveMarket(Mr_plus_2JDtuDJ_Matrix, prefix + "Mat_Mr_plus_2JDtuDJ.mtx");
    Eigen::saveMarket(Inv_Mr_plus_2JDtuDJ_Matrix, prefix + "Mat_Inv_Mr_plus_2JDtuDJ.mtx");
    Eigen::saveMarket(u_Matrix, prefix + "Mat_u.mtx");
    Eigen::saveMarket(uInv_Matrix, prefix + "Mat_uInv.mtx");
    Eigen::saveMarket(uRed_Matrix, prefix + "Mat_uRed.mtx");
    Eigen::saveMarket(uInvRed_Matrix, prefix + "Mat_uInvRed.mtx");
    Eigen::saveMarket(G_Matrix, prefix + "Mat_G.mtx");
    Eigen::saveMarket(Dt_Matrix, prefix + "Mat_Dt.mtx");
    Eigen::saveMarket(JG_Matrix, prefix + "Mat_JG.mtx");
    Eigen::saveMarket(JDt_Matrix, prefix + "Mat_JDt.mtx");
    Eigen::saveMarket(JDtRed_Matrix, prefix + "Mat_JDtRed.mtx");
    Eigen::saveMarketVector(activeRHSVector, prefix + "Vec_activeRHS.mtx");
    Eigen::saveMarketVector(reducedRHSVector, prefix + "Vec_reducedRHS.mtx");
    Eigen::saveMarketVector(pressureRHSVector, prefix + "Vec_pressureRHS.mtx");
    Eigen::saveMarketVector(stressRHSVector, prefix + "Vec_stressRHS.mtx");
}

void
HDK_PolyStokes::Solver::exportMatricesPostSolve(std::string prefix)
{
    Eigen::saveMarketVector(solutionVector, prefix + "solutionVector.mtx");
}

void
HDK_PolyStokes::Solver::exportStats(std::string prefix)
{
    // save a bunch of variables into a matrix
    Eigen::Vector<SolveReal, 27> dimData(
        nCenter,                            // 1
        nFaceX, nFaceY, nFaceZ,             // 2, 3, 4
        nEdgeYZ, nEdgeXZ, nEdgeXY,          // 5, 6, 7
        nActiveVs,                          // 8
        nActiveVx, nActiveVy, nActiveVz,    // 9, 10, 11
        nReducedVs,                         // 12
        nPressures,                         // 13
        nStresses,                          // 14
        nTxx, nTyy, nTzz,                   // 15, 16, 17
        nTyz, nTxz, nTxy,                   // 18, 19, 20
        nTotalDOFs,                         // 21
        nSystemSize,                        // 22
        myThreadCount, myGrainSize,         // 23, 24
        myInteriorRegionCount,              // 25
        dx, dt);                            // 26, 27

    std::cout << prefix << std::endl;
    Eigen::Vector<SolveReal, 6> solveData(
        solveError,             // 1
        solveIterations,        // 2
        solveCPUTime,           // 3
        solveWallclockTime,     // 4
        setupCPUTime,           // 5
        setupWallclockTime);    // 6

    Eigen::saveMarketVector(dimData, prefix + "dimData.mtx");
    Eigen::saveMarketVector(solveData, prefix + "solveData.mtx");
}

void
HDK_PolyStokes::Solver::printStats()
{
    std::cout << "ndofs: " << nTotalDOFs << std::endl;
    std::cout << "Setup CPU time: " << setupCPUTime << " ms" << std::endl;
    std::cout << "Setup Wallclock time: " << setupWallclockTime << " ms" << std::endl;
    std::cout << "Solve CPU time: " << solveCPUTime << " ms" << std::endl;
    std::cout << "Solve Wallclock time: " << solveWallclockTime << " ms" << std::endl;
    std::cout << "Solve Iterations: " << solveIterations << std::endl;
    std::cout << "Solve Error: " << solveError << std::endl;
    std::cout << std::endl;
}

void
HDK_PolyStokes::Solver::readInWarmStart()
{
    Vector inputVector;
    Eigen::loadMarketVector(inputVector, "output_data/test_solutionVector.mtx");
    Eigen::saveMarketVector(inputVector, "output_data/out_test_solutionVector.mtx");

    // the loaded format is velocity - pressure
    exint activeVsOffset = 0;
    exint reducedVsOffset = nActiveVs;
    exint pressureOffset = nActiveVs + nReducedVs;
    exint stressOffset = nActiveVs + nReducedVs + nPressures;

    // copy known values in
    for (exint i = 0; i < nActiveVs; ++i)
        activeGuessVector(i) = inputVector(activeVsOffset + i);
    for (exint i = 0; i < nPressures; ++i)
        pressureGuessVector(i) = inputVector(pressureOffset + i);
    // now fill in the missing dofs
        // todo make sure to include reduced guess vector here
    stressGuessVector = 2. * u_Matrix * (Dt_Matrix.transpose() * activeGuessVector);
}

#define USE_ITERATIVE

HDK_PolyStokes::Solver::SolverResult
HDK_PolyStokes::Solver::solve()
{
    HDK_PolyStokes::Solver::SolverResult solverResult = HDK_PolyStokes::Solver::SolverResult::NOCHANGE;

    if (!checkSolverMatrixCompatibility())
        std::cout << "WARNING: incompatible matrix with solver, ignoring choice of matrix format" << std::endl;

    switch (mySolverType)
    {
        case HDK_PolyStokes_Options::SolverType::EIGEN:
            solverResult = solveEigenCG();
            break;
        case HDK_PolyStokes_Options::SolverType::PCG_DIRECT_PRODUCTS:
            solverResult = solvePCG();
            break;
        case HDK_PolyStokes_Options::SolverType::PCG_MATRIX_VECTOR_PRODUCTS:
            solverResult = solveSPDwithMatrixVectorPCG();
            break;
    }

    return solverResult;
}

bool
HDK_PolyStokes::Solver::checkSolverMatrixCompatibility()
{
    switch (mySolverType)
    {
        case HDK_PolyStokes_Options::SolverType::EIGEN:
            return true;
        case HDK_PolyStokes_Options::SolverType::PCG_DIRECT_PRODUCTS:
            return true;
        case HDK_PolyStokes_Options::SolverType::PCG_MATRIX_VECTOR_PRODUCTS:
            return (myMatrixScheme == HDK_PolyStokes_Options::MatrixScheme::PRESSURE_STRESS);
            break;
        case HDK_PolyStokes_Options::SolverType::SPLIT:
            //return (myMatrixScheme == HDK_PolyStokes_Options::MatrixScheme::PRESSURE_STRESS);
            return true;
            break;
        case HDK_PolyStokes_Options::SolverType::ADMM_FULL:
            //return (myMatrixScheme == HDK_PolyStokes_Options::MatrixScheme::PRESSURE_STRESS);
            return true;
            break;
    }
    return false;
}


HDK_PolyStokes::Solver::SolverResult
HDK_PolyStokes::Solver::solvePCG()
{
    SolverResult result = SolverResult::NOCHANGE;

    std::cout << "myThreadCount: " << myThreadCount << std::endl;
    std::cout << "Eigen threads: " << Eigen::nbThreads() << std::endl;

    // cg temp vars
    Vector tmp_r, tmp_z, tmp_p, tmp_Ap;

    std::clock_t solveCPUClockStart = std::clock();
    std::chrono::time_point<std::chrono::high_resolution_clock> solveWallClockStart = std::chrono::high_resolution_clock::now();

    double totalTimerAapply, totalTimerOther;

    solutionVector.setZero();
    solveIterations = pcg(solutionVector, 
        A, 
        b, 
        tmp_r, tmp_z, tmp_p, tmp_Ap, 
        presolver,
        totalTimerAapply,
        totalTimerOther,
        mySolverTolerance, 
        mySolverMaxIterations);

    std::cout << "timerAapply: " << totalTimerAapply << std::endl << std::flush;
    std::cout << "timerOther: " << totalTimerOther << std::endl << std::flush;

    std::clock_t solveCPUClockEnd = std::clock();
    std::chrono::time_point<std::chrono::high_resolution_clock> solveWallClockEnd = std::chrono::high_resolution_clock::now();

    solveCPUTime = 1000.0 * (solveCPUClockEnd - solveCPUClockStart) / CLOCKS_PER_SEC;
    solveWallclockTime = std::chrono::duration<double, std::milli>(solveWallClockEnd - solveWallClockStart).count();

    return SolverResult::SUCCESS;
}

HDK_PolyStokes::Solver::SolverResult
HDK_PolyStokes::Solver::solveSPDwithMatrixVectorPCG()
{

    //assembleReducedInvertedBlockIncludingBoundaryStresses();

    // first setup our matrix A application
    applyMatrix = new ApplyPressureStressMatrix();
    //applyMatrix->setupDirect(A);
    applyMatrix->setupMatrixVectorProducts(
        dt,
        invDt,
        McInv_Matrix,
        Inv_Mr_plus_2JDtuDJ_Matrix, // should this be BInv (the one with boundary stresses)?
        uInv_Matrix,
        G_Matrix,
        JG_Matrix,
        Dt_Matrix,
        JDt_Matrix
    );

    SolverResult result = SolverResult::NOCHANGE;

    //std::cout << "myThreadCount: " << myThreadCount << std::endl;
    //std::cout << "Eigen threads: " << Eigen::nbThreads() << std::endl;

    // cg temp vars
    Vector tmp_r, tmp_z, tmp_p, tmp_Ap;

    std::clock_t solveCPUClockStart = std::clock();
    std::chrono::time_point<std::chrono::high_resolution_clock> solveWallClockStart = std::chrono::high_resolution_clock::now();
    
    double totalTimerAapply, totalTimerOther;

    solutionVector.setZero();
    solveIterations = pcg_external_matrix_A(solutionVector,
        applyMatrix, 
        b,
        tmp_r, tmp_z, tmp_p, tmp_Ap,
        presolver,
        totalTimerAapply,
        totalTimerOther,
        solveError,
        mySolverTolerance,
        mySolverMaxIterations);

    //std::cout << "timerAapply: " << totalTimerAapply << std::endl << std::flush;
    //std::cout << "timerOther: " << totalTimerOther << std::endl << std::flush;

    // have minres as a backup
    if (solveIterations == mySolverMaxIterations)
    {
        std::cout << "CG failed to converge. Using BiCGStab..." << std::endl;

        solutionVector.setZero();
        solveIterations = bicgstab_external_matrix_A(solutionVector,
            applyMatrix,
            b,
            tmp_r, tmp_z, tmp_p, tmp_Ap,
            presolver,
            totalTimerAapply,
            totalTimerOther,
            solveError,
            mySolverTolerance,
            mySolverMaxIterations);
    }


    std::clock_t solveCPUClockEnd = std::clock();
    std::chrono::time_point<std::chrono::high_resolution_clock> solveWallClockEnd = std::chrono::high_resolution_clock::now();

    solveCPUTime = 1000.0 * (solveCPUClockEnd - solveCPUClockStart) / CLOCKS_PER_SEC;
    solveWallclockTime = std::chrono::duration<double, std::milli>(solveWallClockEnd - solveWallClockStart).count();

    if (solveIterations == mySolverMaxIterations)
        return SolverResult::NOCONVERGE;

    return SolverResult::SUCCESS;
}

HDK_PolyStokes::Solver::SolverResult
HDK_PolyStokes::Solver::solveEigenCG()
{
    SolverResult result = SolverResult::NOCHANGE;

#ifdef USE_ITERATIVE
    std::cout << "Starting solve" << std::endl;
    std::cout << "myThreadCount: " << myThreadCount << std::endl;
    std::cout << "Eigen threads: " << Eigen::nbThreads() << std::endl;
    // lets just use bicgstab for now
    Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper> solver;
    //Eigen::BiCGSTAB<SparseMatrix> solver;
    //Eigen::GMRES<SparseMatrix> solver;
    solver.compute(A);
    solver.setMaxIterations(mySolverMaxIterations);
    solver.setTolerance(mySolverTolerance);

    std::clock_t solveCPUClockStart = std::clock();
    std::chrono::time_point<std::chrono::high_resolution_clock> solveWallClockStart = std::chrono::high_resolution_clock::now();

    solutionVector = solver.solveWithGuess(b, guessVector);

    std::clock_t solveCPUClockEnd = std::clock();
    std::chrono::time_point<std::chrono::high_resolution_clock> solveWallClockEnd = std::chrono::high_resolution_clock::now();

    solveIterations = solver.iterations();
    solveError = solver.error();
    solveCPUTime = 1000.0 * (solveCPUClockEnd - solveCPUClockStart) / CLOCKS_PER_SEC;
    solveWallclockTime = std::chrono::duration<double, std::milli>(solveWallClockEnd - solveWallClockStart).count();

#endif
#ifdef USE_LDLT
    Eigen::SimplicialLDLT<SparseMatrix> solver;
    solver.compute(A);
#endif
#ifdef USE_SPARSELU
    Eigen::SparseLU<SparseMatrix> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
#endif

    if (solver.info() == Eigen::Success)
        return SolverResult::SUCCESS;
    else if (solver.info() == Eigen::NoConvergence)
        return SolverResult::NOCONVERGE;
    // other options: NumericalIssue, InvalidInput

    return SolverResult::NOCONVERGE;
}

void
HDK_PolyStokes::Solver::extractResiduals()
{
    residual = A * solutionVector - b;

    // todo generalize this to any matrix layout
    // right now this is only for pressure-velocity 
    Vector ux_residuals = residual.segment(0, nActiveVx);
    Vector uy_residuals = residual.segment(nActiveVx, nActiveVy);
    Vector uz_residuals = residual.segment(nActiveVx+nActiveVy, nActiveVz);
    Vector p_residuals = residual.segment(nActiveVs+nReducedVs, nPressures);

    writeVectorToField(ux_residuals, velXResidual, SamplingType::FACEX);
    writeVectorToField(uy_residuals, velYResidual, SamplingType::FACEY);
    writeVectorToField(uz_residuals, velZResidual, SamplingType::FACEZ);
    writeVectorToField(p_residuals, pressureResidual, SamplingType::CENTER);
}

void
HDK_PolyStokes::Solver::writeVectorToField(Vector data, SIM_RawField& field, SamplingType type)
{
    using SIM::FieldUtils::faceToCellMap;
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();

    UTparallelForEachNumber(field.field()->numTiles(), [&](const UT_BlockedRange<int>& range)
        {
            UT_VoxelArrayIteratorF vit;
            vit.setConstArray(field.field());

            if (boss->opInterrupt())
                return;

            for (exint i = range.begin(); i != range.end(); ++i)
            {
                vit.myTileStart = i;
                vit.myTileEnd = i + 1;
                vit.rewind();

                for (; !vit.atEnd(); vit.advance())
                {
                    UT_Vector3I coord(vit.x(), vit.y(), vit.z());

                    Index activeIndex = gridLocationToActiveIndex(coord, type);
                    if (activeIndex >= 0)
                    {
                        SolveReal val = data(activeIndex);
                        setFieldValue(field, coord, val);
                    }
                }
            }
        });
}

void
HDK_PolyStokes::Solver::setupClockStart()
{
    setupCPUClockStart = std::clock();
    setupWallClockStart = std::chrono::high_resolution_clock::now();
}

void
HDK_PolyStokes::Solver::setupClockEnd()
{
    setupCPUClockEnd = std::clock();
    setupWallClockEnd = std::chrono::high_resolution_clock::now();

    setupCPUTime = 1000.0 * (setupCPUClockEnd - setupCPUClockStart) / CLOCKS_PER_SEC;
    setupWallclockTime = std::chrono::duration<double, std::milli>(setupWallClockEnd - setupWallClockStart).count();
}

void
HDK_PolyStokes::Solver::applySolutionToVelocity(SIM_RawField& velocity, SIM_RawField& validFaces, const int axis)
{
    using SIM::FieldUtils::faceToCellMap;
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();

    exint validFaceCount = 0;
    // todo these should be functions, used somewhere else
    // create offsets
    exint activeVsOffset = 0;
    exint reducedVsOffset = nActiveVs;
    exint pressureOffset = nActiveVs + nReducedVs;
    exint stressOffset = nActiveVs + nReducedVs + nPressures;

    UTparallelForEachNumber(validFaces.field()->numTiles(), [&](const UT_BlockedRange<int>& range)
        {
            UT_VoxelArrayIteratorF vit;
            vit.setConstArray(validFaces.field());

            if (boss->opInterrupt())
                return;

            for (exint i = range.begin(); i != range.end(); ++i)
            {
                vit.myTileStart = i;
                vit.myTileEnd = i + 1;
                vit.rewind();

                // Tiles for the index grid must have already been uncompressed
                if (!vit.isTileConstant() || vit.getValue() == (float)ValidFlag::VALID_FACE)
                {
                    for (; !vit.atEnd(); vit.advance())
                    {
                        UT_Vector3I face(vit.x(), vit.y(), vit.z());

                        if (vit.getValue() == (float)ValidFlag::INVALID_FACE) {
                            assert(getFieldValue(*liquidWeights[(uint)faceAxisToSamplingType(axis)], face) == 0.);
                            continue;
                        }
                        else {
                            assert(vit.getValue() == (float)ValidFlag::VALID_FACE);
                            assert(getFieldValue(*liquidWeights[(uint)faceAxisToSamplingType(axis)], face) > 0.);
                            validFaceCount++;
                        }

                        exint faceLabel = getFieldValue(*faceLabels(axis), face);
                        exint localActiveFaceIndex = getFieldValue(*faceActiveIndices(axis), face);
                        exint activeFaceIndex = faceVelocityDOF(localActiveFaceIndex, axis);
                        exint reducedFaceIndex = getFieldValue(*faceReducedIndices(axis), face);

                        assert(getFieldValue(*liquidWeights[(uint)faceAxisToSamplingType(axis)], face) > 0.);

                        SolveReal localVelocity = 0.;
                        if (reducedFaceIndex >= 0)
                        { // these should be given by the affine DOFs
                            UT_Vector3T<SolveReal> offset;
                            UT_Vector3T<SolveReal> facePos(face);
                            facePos[axis] -= 0.5;
                            facePos *= dx;
                            offset = facePos - reducedRegionCOM[reducedFaceIndex];

                            RowVector Cx = buildConversionCoefficients(offset, axis);

                            // unchanged value for debugging
                            //ColumnVector reducedVector = reducedRegionBestFitVectors(reducedFaceIndex);
                            //localVelocity = reducedVector.dot(Cx);

                            ColumnVector reducedVector = ColumnVector::Zero();
                            for (int n = 0; n < REDUCED_DOF; ++n)
                            {
                                reducedVector(n) = solutionVector(reducedVsOffset + REDUCED_DOF * reducedFaceIndex + n);
                            }
                            localVelocity = reducedVector.dot(Cx);
                        }
                        else if (localActiveFaceIndex >= 0)
                        { // we solved for velocity on these faces
                            localVelocity = solutionVector(activeVsOffset + activeFaceIndex);
                        }
                        else if (faceLabel == MaterialLabels::SOLID)
                        {
                            localVelocity = getFieldValue(*myCollisionVelocityField->getField(axis), face);
                        }

                        setFieldValue(velocity, face, localVelocity);
                    }
                }
            }
        });
}

void
HDK_PolyStokes::Solver::printAllData()
{
    printIndexData(centerLabels,      SamplingType::CENTER,   "centerLabels");
    printIndexData(faceXLabels,       SamplingType::FACEX,    "faceXLabels");
    printIndexData(faceYLabels,       SamplingType::FACEY,    "faceYLabels");
    printIndexData(faceZLabels,       SamplingType::FACEZ,    "faceZLabels");
    printIndexData(edgeXYLabels,      SamplingType::EDGEXY,   "edgeXYLabels");
    printIndexData(edgeXZLabels,      SamplingType::EDGEXZ,   "edgeXZLabels");
    printIndexData(edgeYZLabels,      SamplingType::EDGEYZ,   "edgeYZLabels");

    printIndexData(centerReducedIndices, SamplingType::CENTER, "centerReducedIndices");
    printIndexData(faceXReducedIndices, SamplingType::FACEX,  "faceXReducedIndices");
    printIndexData(faceYReducedIndices, SamplingType::FACEY,  "faceYReducedIndices");
    printIndexData(faceZReducedIndices, SamplingType::FACEZ,  "faceZReducedIndices");
    printIndexData(edgeXYReducedIndices, SamplingType::EDGEXY, "edgeXYReducedIndices");
    printIndexData(edgeXZReducedIndices, SamplingType::EDGEXZ, "edgeXZReducedIndices");
    printIndexData(edgeYZReducedIndices, SamplingType::EDGEYZ, "edgeYZReducedIndices");

    printIndexData(centerActiveIndices, SamplingType::CENTER, "centerActiveIndices");
    printIndexData(faceXActiveIndices, SamplingType::FACEX, "faceXActiveIndices");
    printIndexData(faceYActiveIndices, SamplingType::FACEY, "faceYActiveIndices");
    printIndexData(faceZActiveIndices, SamplingType::FACEZ, "faceZActiveIndices");
    printIndexData(edgeXYActiveIndices, SamplingType::EDGEXY, "edgeXYActiveIndices");
    printIndexData(edgeXZActiveIndices, SamplingType::EDGEXZ, "edgeXZActiveIndices");
    printIndexData(edgeYZActiveIndices, SamplingType::EDGEYZ, "edgeYZActiveIndices");

    printData(centerFluidWeights,     SamplingType::CENTER,   "centerFluidWeights");
    printData(faceXFluidWeights,      SamplingType::FACEX,    "faceXFluidWeights");
    printData(faceYFluidWeights,      SamplingType::FACEY,    "faceYFluidWeights");
    printData(faceZFluidWeights,      SamplingType::FACEZ,    "faceZFluidWeights");
    printData(edgeXYFluidWeights,     SamplingType::EDGEXY,   "edgeXYFluidWeights");
    printData(edgeXZFluidWeights,     SamplingType::EDGEXZ,   "edgeXZFluidWeights");
    printData(edgeYZFluidWeights,     SamplingType::EDGEYZ,   "edgeYZFluidWeights");

    printData(centerLiquidWeights,    SamplingType::CENTER,   "centerLiquidWeights");
    printData(faceXLiquidWeights,     SamplingType::FACEX,    "faceXLiquidWeights");
    printData(faceYLiquidWeights,     SamplingType::FACEY,    "faceYLiquidWeights");
    printData(faceZLiquidWeights,     SamplingType::FACEZ,    "faceZLiquidWeights");
    printData(edgeXYLiquidWeights,    SamplingType::EDGEXY,   "edgeXYLiquidWeights");
    printData(edgeXZLiquidWeights,    SamplingType::EDGEXZ,   "edgeXZLiquidWeights");
    printData(edgeYZLiquidWeights,    SamplingType::EDGEYZ,   "edgeYZLiquidWeights");

    printArrayData(reducedRegionCOM, "reducedRegionCenterOfMass");
}

void
HDK_PolyStokes::Solver::printResidualData()
{
    printData(pressureResidual, SamplingType::CENTER, "pressureResidual");

    //printData(velXResidual, SamplingType::FACEX, "velXResidual");
    //printData(velYResidual, SamplingType::FACEY, "velYResidual");
    //printData(velZResidual, SamplingType::FACEZ, "velZResidual");
}

// todo  fix this multithreading
void
HDK_PolyStokes::Solver::printData(SIM_RawField& data, SamplingType sampling, const char* name)
{
    UT_Interrupt* boss = UTgetInterrupt();

    SIM_GeometryCopy* geo = myParent.getOrCreateGeometry(myParent.myObj, name);
    SIM_GeometryAutoWriteLock lock(geo, SIM_DATA_ID_PRESERVE);
    GU_Detail* detail = &lock.getGdp();
    detail->clear();
    GA_RWHandleF pscaleHandle(detail, GA_ATTRIB_POINT, "pscale");
    if (!pscaleHandle.isValid())
    {
        detail->addFloatTuple(GA_ATTRIB_POINT, "pscale", 1, GA_Defaults(0));
        pscaleHandle = GA_RWHandleF(detail, GA_ATTRIB_POINT, "pscale");
        pscaleHandle.bumpDataId();
    }
    GA_RWHandleF dataHandle(detail, GA_ATTRIB_POINT, "data");
    if (!dataHandle.isValid())
    {
        detail->addFloatTuple(GA_ATTRIB_POINT, "data", 1, GA_Defaults(-1.));
        dataHandle = GA_RWHandleF(detail, GA_ATTRIB_POINT, "data");
        dataHandle.bumpDataId();
    }

    // Multithread create a list compiling data with point positions
    using DataEntry = std::pair<UT_Vector3, fpreal>;
    using DataArray = std::vector<DataEntry>;
    DataArray dataArray;

    UT_Array<DataArray> parallelDataArray;
    parallelDataArray.setSize(myThreadCount);

    UT_ThreadedAlgorithm compilePointDataAlgo;
    compilePointDataAlgo.run([&](const UT_JobInfo& info)
        {
            DataArray &localDataArray = parallelDataArray[info.job()];
            UT_VoxelArrayIteratorF vit;
            UT_VoxelTileIteratorF vitt;
            vit.setConstArray(data.field());
            vit.splitByTile(info);
            for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
            {
                if (boss->opInterrupt()) break;
                vitt.setTile(vit);
                for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                {
                    UT_Vector3I idx(vitt.x(), vitt.y(), vitt.z());
                    UT_Vector3 pos = orig + ((UT_Vector3)idx + SamplingOffset(sampling)) * dx;
                    localDataArray.emplace_back(pos, vitt.getValue());
                }
            }

            return 0;
        });

    // Save position and data into the detail
    for (UT_Array<DataArray>::iterator it = parallelDataArray.begin(); it != parallelDataArray.end(); it++)
    {
        for (const auto& pointDataPair : *it)
        {
            GA_Offset pointOffset = detail->appendPoint();
            pscaleHandle.set(pointOffset, dx);
            dataHandle.set(pointOffset, pointDataPair.second);
            detail->setPos3(pointOffset, pointDataPair.first);
        }
    }

    pscaleHandle.bumpDataId();
    dataHandle.bumpDataId();
    detail->getAttributes().bumpAllDataIds(GA_ATTRIB_POINT);
}

void
HDK_PolyStokes::Solver::printIndexData(SIM_RawIndexField& data, SamplingType sampling, const char* name)
{
    UT_Interrupt* boss = UTgetInterrupt();

    SIM_GeometryCopy* geo = myParent.getOrCreateGeometry(myParent.myObj, name);
    SIM_GeometryAutoWriteLock lock(geo, SIM_DATA_ID_PRESERVE);
    GU_Detail* detail = &lock.getGdp();
    detail->clear();
    GA_RWHandleF pscaleHandle(detail, GA_ATTRIB_POINT, "pscale");
    if (!pscaleHandle.isValid())
    {
        detail->addFloatTuple(GA_ATTRIB_POINT, "pscale", 1, GA_Defaults(0));
        pscaleHandle = GA_RWHandleF(detail, GA_ATTRIB_POINT, "pscale");
        pscaleHandle.bumpDataId();
    }
    GA_RWHandleF dataHandle(detail, GA_ATTRIB_POINT, "data");
    if (!dataHandle.isValid())
    {
        detail->addFloatTuple(GA_ATTRIB_POINT, "data", 1, GA_Defaults(-1.));
        dataHandle = GA_RWHandleF(detail, GA_ATTRIB_POINT, "data");
        dataHandle.bumpDataId();
    }

    // Multithread create a list compiling data with point positions
    using DataEntry = std::pair<UT_Vector3, int>;
    using DataArray = std::vector<DataEntry>;
    DataArray dataArray;

    UT_Array<DataArray> parallelDataArray;
    parallelDataArray.setSize(myThreadCount);

    UT_ThreadedAlgorithm compilePointDataAlgo;
    compilePointDataAlgo.run([&](const UT_JobInfo& info)
        {
            DataArray& localDataArray = parallelDataArray[info.job()];
            UT_VoxelArrayIteratorI vit;
            UT_VoxelTileIteratorI vitt;
            vit.setConstArray(data.field());
            vit.splitByTile(info);
            for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
            {
                if (boss->opInterrupt()) break;
                vitt.setTile(vit);
                for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                {
                    UT_Vector3I idx(vitt.x(), vitt.y(), vitt.z());
                    UT_Vector3 pos = orig + ((UT_Vector3)idx + SamplingOffset(sampling)) * dx;
                    localDataArray.emplace_back(pos, vitt.getValue());
                }
            }

            return 0;
        });

    // Save position and data into the detail
    for (UT_Array<DataArray>::iterator it = parallelDataArray.begin(); it != parallelDataArray.end(); it++)
    {
        for (const auto& pointDataPair : *it)
        {
            GA_Offset pointOffset = detail->appendPoint();
            pscaleHandle.set(pointOffset, dx);
            dataHandle.set(pointOffset, pointDataPair.second);
            detail->setPos3(pointOffset, pointDataPair.first);
        }
    }

    pscaleHandle.bumpDataId();
    dataHandle.bumpDataId();
    detail->getAttributes().bumpAllDataIds(GA_ATTRIB_POINT);
}

void
HDK_PolyStokes::Solver::printArrayData(UT_Array<UT_Vector3T<SolveReal>>& pos, const char* name)
{
    UT_Interrupt* boss = UTgetInterrupt();

    SIM_GeometryCopy* geo = myParent.getOrCreateGeometry(myParent.myObj, name);
    SIM_GeometryAutoWriteLock lock(geo, SIM_DATA_ID_PRESERVE);
    GU_Detail* detail = &lock.getGdp();
    detail->clear();
    GA_RWHandleF pscaleHandle(detail, GA_ATTRIB_POINT, "pscale");
    if (!pscaleHandle.isValid())
    {
        detail->addFloatTuple(GA_ATTRIB_POINT, "pscale", 1, GA_Defaults(0));
        pscaleHandle = GA_RWHandleF(detail, GA_ATTRIB_POINT, "pscale");
        pscaleHandle.bumpDataId();
    }
    GA_RWHandleF dataHandle(detail, GA_ATTRIB_POINT, "data");
    if (!dataHandle.isValid())
    {
        detail->addFloatTuple(GA_ATTRIB_POINT, "data", 1, GA_Defaults(-1.));
        dataHandle = GA_RWHandleF(detail, GA_ATTRIB_POINT, "data");
        dataHandle.bumpDataId();
    }

    // Save position and data into the detail
    for (UT_Array<UT_Vector3T<SolveReal>>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        GA_Offset pointOffset = detail->appendPoint();
        pscaleHandle.set(pointOffset, dx);
        dataHandle.set(pointOffset, it-pos.begin());
        detail->setPos3(pointOffset, *it + (UT_Vector3T<SolveReal>)orig);
    }

    pscaleHandle.bumpDataId();
    dataHandle.bumpDataId();
    detail->getAttributes().bumpAllDataIds(GA_ATTRIB_POINT);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// called in computeCenterOfMasses()
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
HDK_PolyStokes::Solver::buildReducedRegionCOM(
    UT_Array<UT_Array<UT_Vector3T<SolveReal>>>& parallelInteriorRegionCOM,
    UT_Array<UT_Array<SolveReal>>& parallelInteriorRegionCellCount)
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildInteriorRegionCOMAlgorithm;
    buildInteriorRegionCOMAlgorithm.run([&](const UT_JobInfo& info)
        {
            UT_Array<UT_Vector3T<SolveReal>>& localInteriorRegionCOM = parallelInteriorRegionCOM[info.job()];
            UT_Array<SolveReal>& localInteriorRegionCellCount = parallelInteriorRegionCellCount[info.job()];

            UT_VoxelArrayIteratorI vit;
            vit.setConstArray(centerLabels.field());
            vit.splitByTile(info);

            UT_VoxelTileIteratorI vitt;

            for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
            {
                if (boss->opInterrupt())
                    break;

                if (!vit.isTileConstant() || isReduced(vit.getValue()) )
                {
                    vitt.setTile(vit);

                    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                    {
                        if ( isReduced(vitt.getValue()) )
                        {
                            UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

                            exint interiorRegion = getFieldValue(centerReducedIndices, cell);
                            assert(interiorRegion >= 0);

                            ++localInteriorRegionCellCount[interiorRegion];
                            localInteriorRegionCOM[interiorRegion] += UT_Vector3T<SolveReal>(cell);
                        }
                    }
                }
            }

            return 0;
        });
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// called in computeLeastSquaresFits()
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
HDK_PolyStokes::Solver::buildInteriorBestFitSystems(
    UT_Array<UT_Array<ReducedMatrix>>& parallelBestFitMatrix,
    UT_Array<UT_Array<ColumnVector>>& parallelBestFitRHS)
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildBestFitSystemAlgorithm;
    buildBestFitSystemAlgorithm.run([&](const UT_JobInfo& info)
        {
            UT_VoxelArrayIteratorI vit;
            vit.setConstArray(centerReducedIndices.field());
            vit.splitByTile(info);

            UT_VoxelTileIteratorI vitt;

            UT_Array<ReducedMatrix>& localBestFitMatrix = parallelBestFitMatrix[info.job()];
            UT_Array<ColumnVector>& localBestFitRHS = parallelBestFitRHS[info.job()];

            for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
            {
                if (boss->opInterrupt())
                    break;

                if (!vit.isTileConstant() || vit.getValue() >= 0)
                {
                    vitt.setTile(vit);

                    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                    {
                        UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

                        exint interiorRegion = vitt.getValue();
                        if (interiorRegion >= 0)
                        {
                            for (int axis : {0, 1, 2})
                                for (int direction : {0, 1})
                                {
                                    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);
                                    assert(adjacentCell[axis] >= 0 && adjacentCell[axis] < centerReducedIndices.getVoxelRes()[axis]);

                                    if ( isActive(centerLabels(adjacentCell[0], adjacentCell[1], adjacentCell[2])) )
                                    {
                                        UT_Vector3T<SolveReal> offset(cell);

                                        if (direction == 0)
                                            offset[axis] -= .5;
                                        else
                                            offset[axis] += .5;

                                        offset *= dx;
                                        offset -= reducedRegionCOM[interiorRegion];

                                        ColumnVector columnVector = buildConversionCoefficients(offset, axis);
                                        localBestFitMatrix[interiorRegion] += columnVector * columnVector.transpose();
                                        localBestFitRHS[interiorRegion] += getFieldValue(*myVelocityField->getField(axis), cellToFaceMap(cell, axis, direction)) * columnVector;
                                    }
                                }
                        }
                    }
                }
            }

            return 0;
        });
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// called in computeReducedMassMatrices()
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
HDK_PolyStokes::Solver::buildReducedMassMatrixSystems(
    UT_Array<UT_Array<ReducedMatrix>>& parallelReducedMassMatrix)
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildInteriorMassMatrixSystemAlgorithm;
    buildInteriorMassMatrixSystemAlgorithm.run([&](const UT_JobInfo& info)
        {
            UT_VoxelArrayIteratorI vit;
            vit.setConstArray(centerReducedIndices.field());
            vit.splitByTile(info);

            UT_VoxelTileIteratorI vitt;

            UT_Array<ReducedMatrix>& localMassMatrix = parallelReducedMassMatrix[info.job()];

            for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
            {
                if (boss->opInterrupt())
                    break;

                if (!vit.isTileConstant() || vit.getValue() >= 0)
                {
                    vitt.setTile(vit);

                    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                    {
                        UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

                        exint interiorRegion = vitt.getValue();
                        if (interiorRegion >= 0)
                        {
                            for (int axis : {0, 1, 2}) {
                                for (int direction : {0, 1})
                                {
                                    bool doApplyFace = false;
                                    if (direction == 0)
                                        doApplyFace = true;
                                    else
                                    {
                                        UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);
                                        assert(adjacentCell[axis] >= 0 && adjacentCell[axis] < centerReducedIndices.getVoxelRes()[axis]);

                                        if ( isActive(getFieldValue(centerLabels, adjacentCell)) )
                                            doApplyFace = true;
                                    }

                                    if (doApplyFace)
                                    {
                                        UT_Vector3T<SolveReal> offset(cell);
                                        UT_Vector3I face = cellToFaceMap(cell, axis, direction);

                                        if (direction == 0)
                                            offset[axis] -= .5;
                                        else
                                            offset[axis] += .5;

                                        offset *= dx;
                                        offset -= reducedRegionCOM[interiorRegion];

                                        ColumnVector columnVector = buildConversionCoefficients(offset, axis);
                                        localMassMatrix[interiorRegion] += myConstantDensity * columnVector * columnVector.transpose();
                                    }
                                }
                            }
                        }
                    }
                }
            }

            return 0;
        });
}

void
HDK_PolyStokes::Solver::buildReducedViscosityMatrixSystemsInteriorOnly(
    UT_Array<UT_Array<ReducedMatrix>>& parallelReducedViscosityMatrix)
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::faceToCellMap;
    using SIM::FieldUtils::faceToEdgeMap;
    using SIM::FieldUtils::edgeToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    for (int faceAxis : {0, 1, 2})
    {
        UT_Interrupt* boss = UTgetInterrupt();
        UT_ThreadedAlgorithm buildInteriorViscositySystemAlgorithm;
        buildInteriorViscositySystemAlgorithm.run([&](const UT_JobInfo& info)
            {
                UT_Array<ReducedMatrix>& localInteriorViscosityMatrix = parallelReducedViscosityMatrix[info.job()];

                UT_VoxelArrayIteratorI vit;
                SIM_RawIndexField tileFaceIndex;
                switch (faceAxis)
                {
                case 0:
                    tileFaceIndex = faceXReducedIndices;
                    break;
                case 1:
                    tileFaceIndex = faceYReducedIndices;
                    break;
                case 2:
                    tileFaceIndex = faceZReducedIndices;
                    break;
                }
                vit.setConstArray(tileFaceIndex.field());
                vit.splitByTile(info);

                UT_VoxelTileIteratorI vitt;
                for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
                {
                    if (boss->opInterrupt())
                        break;

                    if (!vit.isTileConstant() || vit.getValue() >= 0)
                    {
                        vitt.setTile(vit);
                        for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                        {
                            UT_Vector3I face(vitt.x(), vitt.y(), vitt.z());

                            exint selfInteriorIndex = vitt.getValue();
                            if (selfInteriorIndex < 0) continue;
                            assert(selfInteriorIndex < myInteriorRegionCount);

                            // build cell-centered stress terms
                            for (int divergenceDirection : {0, 1})
                            {
                                UT_Vector3I cell = faceToCellMap(face, faceAxis, divergenceDirection);
                                // first ignore stencils of anything that isn't an interior dof
                                if ( !isReduced(centerLabels(cell.x(), cell.y(), cell.z())) )
                                    continue;
                                // oob check
                                if (cell[faceAxis] < 0 || cell[faceAxis] >= tileFaceIndex.getVoxelRes()[faceAxis])
                                    continue;
                                SolveReal divergenceSign = (divergenceDirection == 0) ? -1 : 1;
                                // grab viscosity
                                UT_Vector3 point;
                                centerLabels.indexToPos(cell[0], cell[1], cell[2], point);
                                SolveReal visc = getLocalViscosity(point);

                                for (int gradientDirection : {0, 1})
                                {
                                    UT_Vector3I adjacentFace = cellToFaceMap(cell, faceAxis, gradientDirection);
                                    SolveReal gradientSign = (gradientDirection == 0) ? -1. : 1.;
                                    // we don't need the fluid and liquid volumes since we assume the interior regions
                                    // are sufficiently far from the surface
                                    SolveReal contribution = -1. * divergenceSign
                                        * gradientSign
                                        * visc / (dx * dx);
                                    SIM_RawIndexField* adjacentTileFaceIndex;
                                    switch (faceAxis)
                                    {
                                    case 0:
                                        adjacentTileFaceIndex = &faceXReducedIndices;
                                        break;
                                    case 1:
                                        adjacentTileFaceIndex = &faceYReducedIndices;
                                        break;
                                    case 2:
                                        adjacentTileFaceIndex = &faceZReducedIndices;
                                        break;
                                    }
                                    exint adjacentInteriorIndex = getFieldValue(*adjacentTileFaceIndex, adjacentFace);
                                    if (adjacentInteriorIndex < 0) continue;
                                    //assert(adjacentInteriorIndex < myInteriorRegionCount);

                                    
                                    // compute offsets
                                    UT_Vector3T<SolveReal> selfOffset;
                                    UT_Vector3T<SolveReal> adjOffset;

                                    {
                                        UT_Vector3T<SolveReal> facePos(face);
                                        facePos[faceAxis] -= 0.5;
                                        facePos *= dx;
                                        selfOffset = facePos - reducedRegionCOM[selfInteriorIndex];
                                    }
                                    {
                                        UT_Vector3T<SolveReal> facePos(adjacentFace);
                                        facePos[faceAxis] -= 0.5;
                                        facePos *= dx;
                                        adjOffset = facePos - reducedRegionCOM[adjacentInteriorIndex];
                                    }

                                    RowVector rowVec = buildConversionCoefficients(adjOffset, faceAxis).transpose();
                                    ColumnVector colVec = buildConversionCoefficients(selfOffset, faceAxis);

                                    assert(adjacentInteriorIndex == selfInteriorIndex);
                                    localInteriorViscosityMatrix[selfInteriorIndex] += contribution * colVec * rowVec;
                                }
                            }
                            
                            // build edge-centered stress terms
                            for (int edgeAxis : {0, 1, 2})
                            {
                                if (edgeAxis == faceAxis) continue;
                                for (int divergenceDirection : {0, 1})
                                {
                                    SolveReal divergenceSign = (divergenceDirection == 0) ? -1 : 1;
                                    UT_Vector3I edge = faceToEdgeMap(face, faceAxis, edgeAxis, divergenceDirection);

                                    // grab viscosity
                                    UT_Vector3 point;
                                    edgeLabels(edgeAxis)->indexToPos(edge[0], edge[1], edge[2], point);
                                    float visc = getLocalViscosity(point);

                                    bool edgeIsInterior = false;
                                    edgeIsInterior = isReducedButNotBoundary(getFieldValue(*edgeLabels(edgeAxis), edge));

                                    // first ignore stencils which isn't an interior dof
                                    if (!edgeIsInterior)
                                        continue;

                                    for (int gradientAxis : {0, 1, 2})
                                    {
                                        if (gradientAxis == edgeAxis) continue;
                                        int adjacentFaceAxis = 3 - gradientAxis - edgeAxis;
                                        SIM_RawIndexField* adjacentTileFaceIndex;
                                        switch (adjacentFaceAxis)
                                        {
                                        case 0:
                                            adjacentTileFaceIndex = &faceXReducedIndices;
                                            break;
                                        case 1:
                                            adjacentTileFaceIndex = &faceYReducedIndices;
                                            break;
                                        case 2:
                                            adjacentTileFaceIndex = &faceZReducedIndices;
                                            break;
                                        }

                                        for (int gradientDirection : {0, 1})
                                        {
                                            UT_Vector3I adjacentFace = edgeToFaceMap(edge, edgeAxis, adjacentFaceAxis, gradientDirection);
                                            SolveReal gradientSign = (gradientDirection == 0) ? -1 : 1;
                                            SolveReal contribution = -0.5 * divergenceSign
                                                * gradientSign
                                                * visc / (dx * dx);

                                            // TOCHECK make sure this boundary check is right
                                            //if (adjacentFace[gradientAxis] < 0 || adjacentFace[gradientAxis] >= activeFaceIndices[faceAxis].getVoxelRes()[gradientAxis])
                                            //    continue;
                                            int adjacentInteriorIndex = getFieldValue(*adjacentTileFaceIndex, adjacentFace);
                                            if (adjacentInteriorIndex < 0) continue;

                                            // compute offsets
                                            UT_Vector3T<SolveReal> selfOffset;
                                            UT_Vector3T<SolveReal> adjOffset;

                                            {
                                                UT_Vector3T<SolveReal> facePos(face);
                                                facePos[faceAxis] -= 0.5;
                                                facePos *= dx;
                                                selfOffset = facePos - reducedRegionCOM[selfInteriorIndex];
                                            }
                                            {
                                                UT_Vector3T<SolveReal> facePos(adjacentFace);
                                                facePos[adjacentFaceAxis] -= 0.5;
                                                facePos *= dx;
                                                adjOffset = facePos - reducedRegionCOM[adjacentInteriorIndex];
                                            }

                                            RowVector rowVec = buildConversionCoefficients(adjOffset, adjacentFaceAxis).transpose();
                                            ColumnVector colVec = buildConversionCoefficients(selfOffset, faceAxis);

                                            //assert(adjacentInteriorIndex == selfInteriorIndex);
                                            localInteriorViscosityMatrix[selfInteriorIndex] += contribution * colVec * rowVec;
                                        }
                                    }
                                }
                            }
                            
                        }
                    }
                }

                return 0;
            });
    }

    return;
}

void
HDK_PolyStokes::Solver::buildReducedViscosityMatrixSystems(
    UT_Array<UT_Array<ReducedMatrix>>& parallelReducedViscosityMatrix)
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::faceToCellMap;
    using SIM::FieldUtils::faceToEdgeMap;
    using SIM::FieldUtils::edgeToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    for (int faceAxis : {0, 1, 2})
    {
        UT_Interrupt* boss = UTgetInterrupt();
        UT_ThreadedAlgorithm buildInteriorViscositySystemAlgorithm;
        buildInteriorViscositySystemAlgorithm.run([&](const UT_JobInfo& info)
            {
                UT_Array<ReducedMatrix>& localInteriorViscosityMatrix = parallelReducedViscosityMatrix[info.job()];

                UT_VoxelArrayIteratorI vit;
                SIM_RawIndexField tileFaceIndex;
                switch (faceAxis)
                {
                case 0:
                    tileFaceIndex = faceXReducedIndices;
                    break;
                case 1:
                    tileFaceIndex = faceYReducedIndices;
                    break;
                case 2:
                    tileFaceIndex = faceZReducedIndices;
                    break;
                }
                vit.setConstArray(tileFaceIndex.field());
                vit.splitByTile(info);

                UT_VoxelTileIteratorI vitt;
                for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
                {
                    if (boss->opInterrupt())
                        break;

                    if (!vit.isTileConstant() || vit.getValue() >= 0)
                    {
                        vitt.setTile(vit);
                        for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                        {
                            UT_Vector3I face(vitt.x(), vitt.y(), vitt.z());

                            exint selfInteriorIndex = vitt.getValue();
                            assert(selfInteriorIndex < myInteriorRegionCount);
                            if (selfInteriorIndex >= 0)
                            {
                                // build cell-centered stress terms
                                for (int divergenceDirection : {0, 1})
                                {
                                    UT_Vector3I cell = faceToCellMap(face, faceAxis, divergenceDirection);
                                    // first ignore stencils of anything that isn't an interior dof
                                    //if (centerLabels(cell.x(), cell.y(), cell.z()) != MaterialLabels::REDUCED)
                                    //    continue;
                                    // oob check
                                    if (cell[faceAxis] < 0 || cell[faceAxis] >= tileFaceIndex.getVoxelRes()[faceAxis])
                                        continue;
                                    SolveReal divergenceSign = (divergenceDirection == 0) ? -1 : 1;
                                    // grab viscosity
                                    UT_Vector3 point;
                                    centerLabels.indexToPos(cell[0], cell[1], cell[2], point);
                                    float visc = getLocalViscosity(point);

                                    for (int gradientDirection : {0, 1})
                                    {
                                        UT_Vector3I adjacentFace = cellToFaceMap(cell, faceAxis, gradientDirection);
                                        SolveReal gradientSign = (gradientDirection == 0) ? -1. : 1.;
                                        // we don't need the fluid and liquid volumes since we assume the interior regions
                                        // are sufficiently far from the surface
                                        SolveReal contribution = -1. * divergenceSign
                                            * gradientSign
                                            * visc / (dx * dx);
                                        SIM_RawIndexField adjacentTileFaceIndex;
                                        switch (faceAxis)
                                        {
                                        case 0:
                                            adjacentTileFaceIndex = faceXReducedIndices;
                                            break;
                                        case 1:
                                            adjacentTileFaceIndex = faceYReducedIndices;
                                            break;
                                        case 2:
                                            adjacentTileFaceIndex = faceZReducedIndices;
                                            break;
                                        }
                                        exint adjacentInteriorIndex = getFieldValue(adjacentTileFaceIndex, adjacentFace);
                                        assert(adjacentInteriorIndex < myInteriorRegionCount);

                                        // compute offsets
                                        UT_Vector3T<SolveReal> selfOffset;
                                        UT_Vector3T<SolveReal> adjOffset;
                                        if (selfInteriorIndex >= 0)
                                        {
                                            UT_Vector3T<SolveReal> facePos(face);
                                            facePos[faceAxis] -= 0.5;
                                            facePos *= dx;
                                            selfOffset = facePos - reducedRegionCOM[selfInteriorIndex];
                                        }
                                        if (adjacentInteriorIndex >= 0)
                                        {
                                            UT_Vector3T<SolveReal> facePos(adjacentFace);
                                            facePos[faceAxis] -= 0.5;
                                            facePos *= dx;
                                            adjOffset = facePos - reducedRegionCOM[adjacentInteriorIndex];
                                        }

                                        RowVector rowVec = buildConversionCoefficients(adjOffset, faceAxis).transpose();
                                        ColumnVector colVec = buildConversionCoefficients(selfOffset, faceAxis);

                                        if (selfInteriorIndex >= 0 && adjacentInteriorIndex >= 0)
                                        {
                                            assert(adjacentInteriorIndex == selfInteriorIndex);
                                            localInteriorViscosityMatrix[selfInteriorIndex] += contribution * colVec * rowVec;
                                        }
                                    }
                                }

                                // build edge-centered stress terms
                                for (int edgeAxis : {0, 1, 2})
                                {
                                    if (edgeAxis == faceAxis) continue;
                                    for (int divergenceDirection : {0, 1})
                                    {
                                        SolveReal divergenceSign = (divergenceDirection == 0) ? -1 : 1;
                                        UT_Vector3I edge = faceToEdgeMap(face, faceAxis, edgeAxis, divergenceDirection);

                                        // grab viscosity
                                        UT_Vector3 point;
                                        edgeLabels(edgeAxis)->indexToPos(edge[0], edge[1], edge[2], point);
                                        float visc = getLocalViscosity(point);

                                        bool edgeIsInterior = false;
                                        // first ignore stencils which isn't an interior dof
                                        //if (!edgeIsInterior)
                                          //  continue;

                                        for (int gradientAxis : {0, 1, 2})
                                        {
                                            if (gradientAxis == edgeAxis) continue;
                                            int adjacentFaceAxis = 3 - gradientAxis - edgeAxis;
                                            SIM_RawIndexField adjacentTileFaceIndex;
                                            switch (adjacentFaceAxis)
                                            {
                                            case 0:
                                                adjacentTileFaceIndex = faceXReducedIndices;
                                                break;
                                            case 1:
                                                adjacentTileFaceIndex = faceYReducedIndices;
                                                break;
                                            case 2:
                                                adjacentTileFaceIndex = faceZReducedIndices;
                                                break;
                                            }

                                            for (int gradientDirection : {0, 1})
                                            {
                                                UT_Vector3I adjacentFace = edgeToFaceMap(edge, edgeAxis, adjacentFaceAxis, gradientDirection);
                                                SolveReal gradientSign = (gradientDirection == 0) ? -1 : 1;
                                                SolveReal contribution = -0.5 * divergenceSign
                                                    * gradientSign
                                                    * visc / (dx * dx);

                                                // TOCHECK make sure this boundary check is right
                                                //if (adjacentFace[gradientAxis] < 0 || adjacentFace[gradientAxis] >= activeFaceIndices[faceAxis].getVoxelRes()[gradientAxis])
                                                //    continue;
                                                int adjacentInteriorIndex = getFieldValue(adjacentTileFaceIndex, adjacentFace);

                                                // compute offsets
                                                UT_Vector3T<SolveReal> selfOffset;
                                                UT_Vector3T<SolveReal> adjOffset;
                                                if (selfInteriorIndex >= 0)
                                                {
                                                    UT_Vector3T<SolveReal> facePos(face);
                                                    facePos[faceAxis] -= 0.5;
                                                    facePos *= dx;
                                                    selfOffset = facePos - reducedRegionCOM[selfInteriorIndex];
                                                }
                                                if (adjacentInteriorIndex >= 0)
                                                {
                                                    UT_Vector3T<SolveReal> facePos(adjacentFace);
                                                    facePos[adjacentFaceAxis] -= 0.5;
                                                    facePos *= dx;
                                                    adjOffset = facePos - reducedRegionCOM[adjacentInteriorIndex];
                                                }

                                                RowVector rowVec = buildConversionCoefficients(adjOffset, adjacentFaceAxis).transpose();
                                                ColumnVector colVec = buildConversionCoefficients(selfOffset, faceAxis);

                                                if (selfInteriorIndex >= 0 && adjacentInteriorIndex >= 0)
                                                {
                                                    assert(adjacentInteriorIndex == selfInteriorIndex);
                                                    localInteriorViscosityMatrix[selfInteriorIndex] += contribution * colVec * rowVec;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                return 0;
            });
    }

    return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Utils
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
fpreal
HDK_PolyStokes::Solver::getLocalDensity(UT_Vector3I face)
{
    return myConstantDensity;
}

fpreal
HDK_PolyStokes::Solver::getLocalViscosity(UT_Vector3 point)
{
    return myViscosityFieldData->getValue(point);
}

void
HDK_PolyStokes::Solver::copyMaterialLabel(
    SIM_RawIndexField& source,
    SIM_RawIndexField& dest,
    const exint searchValue)
{
    UT_Interrupt* boss = UTgetInterrupt();
    UTparallelForEachNumber(source.field()->numTiles(), [&](const UT_BlockedRange<int>& range)
        {
            UT_VoxelArrayIteratorI vit(source.fieldNC());
            UT_VoxelArrayIteratorI destit(dest.fieldNC());

            if (boss->opInterrupt())
                return;

            for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
            {
                vit.myTileStart = tileNumber;
                vit.myTileEnd = tileNumber + 1;
                vit.rewind();

                destit.myTileStart = tileNumber;
                destit.myTileEnd = tileNumber + 1;
                destit.rewind();

                if (vit.isTileConstant() && vit.getValue() == searchValue)
                {
                    UT_VoxelTile<exint>* tile = destit.getTile();
                    tile->makeConstant(searchValue);
                }
                else
                {
                    for (; !vit.atEnd(); vit.advance())
                    {
                        if (vit.getValue() == searchValue)
                            destit.setValue(searchValue);
                        destit.advance();
                    }
                }
            }
        });
}

void
HDK_PolyStokes::Solver::overwriteIndices(
    SIM_RawIndexField& indices,
    const exint searchValue,
    const exint replaceValue)
{
    UT_Interrupt* boss = UTgetInterrupt();
    UTparallelForEachNumber(indices.field()->numTiles(), [&](const UT_BlockedRange<int>& range)
        {
            UT_VoxelArrayIteratorI vit(indices.fieldNC());

            if (boss->opInterrupt())
                return;

            for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
            {
                vit.myTileStart = tileNumber;
                vit.myTileEnd = tileNumber + 1;
                vit.rewind();

                if (vit.isTileConstant() && vit.getValue() == searchValue)
                {
                    UT_VoxelTile<exint>* tile = vit.getTile();
                    tile->makeConstant(replaceValue);
                }
                else
                {
                    for (; !vit.atEnd(); vit.advance())
                    {
                        if (vit.getValue() == searchValue)
                            vit.setValue(replaceValue);
                    }
                }
            }
        });
}

template<typename Grid>
void
HDK_PolyStokes::Solver::initField(
    Grid& field,
    const SIM_RawField& surface,
    const SIM_FieldSample sample)
{
    UT_Vector3 size = surface.getSize();
    UT_Vector3 orig = surface.getOrig();
    UT_Vector3i res;
    surface.getVoxelRes(res[0], res[1], res[2]);

    field.init(sample, orig, size, res[0], res[1], res[2]);
    field.makeConstant(0.);
}

void
HDK_PolyStokes::Solver::setActiveLayerCells(const UT_Array<UT_Vector3I>& activeCellLayer)
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();

    // Assign active cells
    UTparallelForLightItems(UT_BlockedRange<exint>(0, activeCellLayer.size()), [&](const UT_BlockedRange<exint>& range)
        {
            if (boss->opInterrupt())
                return;

            exint startIndex = range.begin();

            if (startIndex > 0)
            {
                while (startIndex != range.end() && activeCellLayer[startIndex] == activeCellLayer[startIndex - 1])
                    ++startIndex;
            }

            UT_Vector3I oldCell(-1, -1, -1);

            for (exint cellIndex = startIndex; cellIndex != range.end(); ++cellIndex)
            {
                UT_Vector3I cell = activeCellLayer[cellIndex];

                if (cell == oldCell)
                    continue;

                oldCell = cell;

                assert(getFieldValue(centerLabels, cell) == MaterialLabels::GENERICFLUID);
                setFieldValue(centerLabels, cell, MaterialLabels::ACTIVEFLUID);
            }
        });

}

void
HDK_PolyStokes::Solver::findOccupiedIndexTiles(UT_Array<bool>& isTileOccupiedList,
    const UT_Array<UT_Vector3I>& indexCellList,
    const SIM_RawIndexField& indexCellLabels)
{
    UT_Interrupt* boss = UTgetInterrupt();

    const exint tileSize = isTileOccupiedList.entries();

    UT_Array<bool> localIsTileOccupiedList;
    localIsTileOccupiedList.setSize(tileSize);
    localIsTileOccupiedList.constant(false);

    if (boss->opInterrupt())
        return;

    const exint elementSize = indexCellList.entries();
    UTparallelFor(UT_BlockedRange<exint>(0, elementSize), [&](const UT_BlockedRange<exint>& range)
        {
            for (exint i = range.begin(); i != range.end(); ++i)
            {
                if (!(i & 127))
                {
                    if (boss->opInterrupt())
                        break;
                }

                UT_Vector3I cell = indexCellList[i];

                int tileNumber = indexCellLabels.field()->indexToLinearTile(cell[0], cell[1], cell[2]);

                localIsTileOccupiedList[tileNumber] = true;
            }
        });

    for (exint tileNumber = 0; tileNumber < tileSize; ++tileNumber)
    {
        if (localIsTileOccupiedList[tileNumber])
            isTileOccupiedList[tileNumber] = true;
    }
}


#ifdef QUADRATIC_REGIONS

ColumnVector
HDK_PolyStokes::Solver::buildConversionCoefficients(UT_Vector3T<SolveReal> offset, int axis)
{
    ColumnVector vals = ColumnVector::Zero();

    switch (axis) {
    case 0:
        vals << 1., 0., 0.,
            offset[0], offset[1], offset[2],
            offset[0] * offset[0], offset[0] * offset[1], offset[0] * offset[2],
            offset[1] * offset[1], offset[1] * offset[2], offset[2] * offset[2],
            0., 0., 0.,
            0., 0., 0.,
            0., 0., 0.,
            0., 0., 0.,
            0., 0.;
        break;
    case 1:
        vals << 0., 1., 0.,
            0., 0., 0.,
            0., 0., 0.,
            0., 0., 0.,
            offset[0], offset[1], offset[2],
            offset[0] * offset[0], offset[0] * offset[1], offset[0] * offset[2],
            offset[1] * offset[1], offset[1] * offset[2], offset[2] * offset[2],
            0., 0., 0.,
            0., 0.;
        break;
    case 2:
        vals << 0., 0., 1.,
            -offset[2], 0., 0.,
            -2. * offset[0] * offset[2], -1. * offset[1] * offset[2], -0.5 * offset[2] * offset[2],
            0., 0., 0.,
            0., -offset[2], 0.,
            0., -1. * offset[0] * offset[2], 0.,
            -2. * offset[1] * offset[2], -0.5 * offset[2] * offset[2], 0.,
            offset[0], offset[1], offset[0] * offset[0],
            offset[0] * offset[1], offset[1] * offset[1];
        break;
    }

    return vals;
}

#endif

#ifdef AFFINE_REGIONS

ColumnVector
HDK_PolyStokes::Solver::buildConversionCoefficients(UT_Vector3T<SolveReal> offset, int axis)
{
    ColumnVector vals = ColumnVector::Zero();

    switch (axis) {
    case 0:
        vals << 1., 0., 0.,
            offset[0], offset[1], offset[2],
            0., 0., 0.,
            0., 0.;
        break;
    case 1:
        vals << 0., 1., 0.,
            0., 0., 0.,
            offset[0], offset[1], offset[2],
            0., 0.;
        break;
    case 2:
        vals << 0., 0., 1.,
            -offset[2], 0., 0.,
            0., -offset[2], 0.,
            offset[0], offset[1];
        break;
    }

    return vals;
}

#endif