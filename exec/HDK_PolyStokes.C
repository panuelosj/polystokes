/*
 * 
 */
#include "HDK_PolyStokes.h"
#include "HDK_PolyStokesSolver.h"

// tools for hooking into Houdini as a DOP
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <PRM/PRM_Include.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_DopDescription.h>
#include <UT/UT_PerfMonAutoEvent.h>
// gas solver tools
#include <SIM/SIM_FieldUtils.h>
#include <SIM/SIM_Object.h>
#include <GAS/GAS_SubSolver.h>
// for printing out data
#include <SIM/SIM_GeometryCopy.h>
// multithreading
#include <UT/UT_ThreadedAlgorithm.h>
#include <tbb/tbb.h>
//
#include <string>

// Houdini hook
void
initializeSIM(void*)
{
    IMPLEMENT_DATAFACTORY(HDK_PolyStokes);
}

// Constructor and destructor
HDK_PolyStokes::HDK_PolyStokes(const SIM_DataFactory* factory)
    : BaseClass(factory)
{
    setCenterLabels("centerLabels");
    setFaceXLabelsName("faceXLabels");
    setFaceYLabelsName("faceYLabels");
    setFaceZLabelsName("faceZLabels");
    setEdgeXYLabelsName("edgeXYLabels");
    setEdgeXZLabelsName("edgeXZLabels");
    setEdgeYZLabelsName("edgeYZLabels");

    setCenterReducedIndices("centerReducedIndices");
    setFaceXReducedIndicesName("faceXReducedIndices");
    setFaceYReducedIndicesName("faceYReducedIndices");
    setFaceZReducedIndicesName("faceZReducedIndices");
    setEdgeXYReducedIndicesName("edgeXYReducedIndices");
    setEdgeXZReducedIndicesName("edgeXZReducedIndices");
    setEdgeYZReducedIndicesName("edgeYZReducedIndices");

    setCenterActiveIndices("centerActiveIndices");
    setFaceXActiveIndicesName("faceXActiveIndices");
    setFaceYActiveIndicesName("faceYActiveIndices");
    setFaceZActiveIndicesName("faceZActiveIndices");
    setEdgeXYActiveIndicesName("edgeXYActiveIndices");
    setEdgeXZActiveIndicesName("edgeXZActiveIndices");
    setEdgeYZActiveIndicesName("edgeYZActiveIndices");

    setCenterLiquidWeightsName("centerLiquidWeights");
    setFaceXLiquidWeightsName("faceXLiquidWeights");
    setFaceYLiquidWeightsName("faceYLiquidWeights");
    setFaceZLiquidWeightsName("faceZLiquidWeights");
    setEdgeXYLiquidWeightsName("edgeXYLiquidWeights");
    setEdgeXZLiquidWeightsName("edgeXZLiquidWeights");
    setEdgeYZLiquidWeightsName("edgeYZLiquidWeights");

    setCenterFluidWeightsName("centerFluidWeights");
    setFaceXFluidWeightsName("faceXFluidWeights");
    setFaceYFluidWeightsName("faceYFluidWeights");
    setFaceZFluidWeightsName("faceZFluidWeights");
    setEdgeXYFluidWeightsName("edgeXYFluidWeights");
    setEdgeXZFluidWeightsName("edgeXZFluidWeights");
    setEdgeYZFluidWeightsName("edgeYZFluidWeights");

    setReducedRegionCenterOfMassName("reducedRegionCenterOfMass");
}

HDK_PolyStokes::~HDK_PolyStokes()
{
}

// Create HDK Node
const SIM_DopDescription *
HDK_PolyStokes::getDopDescription()
{
    static PRM_Name     theVelocityName(GAS_NAME_VELOCITY, "Velocity Field");
    static PRM_Default  theVelocityNameDefault(0, "vel");

    static PRM_Name     theValidName("valid", "Valid Field");
    static PRM_Default  theValidNameDefault(0, "__valid");

    static PRM_Name     theViscosityName("viscosity", "Viscosity Field");
    static PRM_Default  theViscosityNameDefault(0, "viscosity");

    static PRM_Name     theDensityName(GAS_NAME_DENSITY, "Liquid Density Field");
    static PRM_Default  theDensityNameDefault(0, "massdensity");
    static PRM_Name     theMinDensityName("mindensity", "Min Density");
    static PRM_Name     theMaxDensityName("maxdensity", "Max Density");
    static PRM_Default  theMaxDensityDefault(100000);

    static PRM_Name     thePressureName(GAS_NAME_PRESSURE, "Pressure Field");
    static PRM_Default  thePressureNameDefault(0, "pressure");

    static PRM_Name     theSurfaceName(GAS_NAME_SURFACE, "Surface Field");
    static PRM_Default  theSurfaceNameDefault(0, "surface");

    static PRM_Name     theCollisionName(GAS_NAME_COLLISION, "Solid Collision Field");
    static PRM_Default  theCollisionNameDefault(0, "collision");

    static PRM_Name     theCollisionVelocityName(GAS_NAME_COLLISIONVELOCITY, "Solid Collision Velocity Field");
    static PRM_Default  theCollisionVelocityNameDefault(0, "collisionvel");

    static PRM_Name     theUseInputSurfaceWeightsName("useInputSurfaceWeights", "Use Input Surface Weights");
    static PRM_Name     theSurfaceWeightsName("surfaceweights", "Surface Weights Field");
    static PRM_Default  theSurfaceWeightsNameDefault(0, "surfaceweights");

    static PRM_Name     theUseInputCollisionWeights("useInputCollisionWeights", "Use Input Collision Weights");
    static PRM_Name     theCollisionWeightsName("collisionweights", "Collision Weights Field");
    static PRM_Default  theCollisionWeightsNameDefault(0, "collisionweights");

    static PRM_Name     theActiveLiquidBoundaryLayerSizeName("activeLiquidBoundaryLayerSize", "Active Liquid Boundary Layer Size");
    static PRM_Default  theActiveLiquidBoundaryLayerSizeDefault(2);
    static PRM_Name     theActiveSolidBoundaryLayerSizeName("activeSolidBoundaryLayerSize", "Active Solid Boundary Layer Size");
    static PRM_Default  theActiveSolidBoundaryLayerSizeDefault(2);
    static PRM_Name     theDoReducedRegionsName("doReducedRegions", "Do Reduced Regions");
    static PRM_Name     theDoTileName("doTile", "Do Tile");
    static PRM_Name     theTileSizeName("tileSize", "Reduced Tile Size");
    static PRM_Default  theTileSizeDefault(16);
    static PRM_Name     theTilePaddingName("tilePadding", "Reduced Tile Padding");
    static PRM_Default  theTilePaddingDefault(2);

    static PRM_Name     theExportDataPrefixName("exportDataPrefix", "Export Data Prefix");
    static PRM_Default  theExportDataPrefixDefault(0, "output_data/`opname(\"../..\")`.$FF.");

    static PRM_Name     theDoSolveName("doSolve", "Do Solve");
    static PRM_Name     theKeepNonConvergedResultsName("keepNonConvergedResults", "Keep Non-Converged Results");
    static PRM_Name     theExportMatricesName("exportMatrices", "Export Matrices");
    static PRM_Name     theExportComponentMatricesName("exportComponentMatrices", "Export Component Matrices");
    static PRM_Name     theExportStatsName("exportStats", "Export Stats");
    static PRM_Name     theUseWarmStartName("useWarmStart", "Use Warm Start");

    static PRM_Name 	theToleranceName(SIM_NAME_TOLERANCE, "Solver Tolerance");
    static PRM_Default 	theToleranceDefault(1e-3);

    static PRM_Name 	theMaxSolverIterations("maxSolverIterations", "Max Solver Iterations");
    static PRM_Default 	theMaxSolverIterationsDefault(5000);

    static PRM_Name     theMatrixSchemeName("matrixSetup", "Matrix Setup");
    static PRM_Name     theMatrixSchemeChoices[] =
    {
        PRM_Name("pressurestress", "Pressure Stress SPD Form"),
        //PRM_Name("fulldofs", "All DOFs (Velocity, Pressure, Stress)"),
        //PRM_Name("pressurevelocity", "Pressure Velocity Form"),
        PRM_Name(0)
    };
    static PRM_ChoiceList theMatrixSchemeMenu(PRM_CHOICELIST_SINGLE, theMatrixSchemeChoices);

    static PRM_Name     theSolverTypeName("solverType", "Solver Type");
    static PRM_Name     theSolverTypeChoices[] =
    {
        PRM_Name("pcg_matrix_vector_products", "Preconditioned CG - Factored Matrix Vector Products"),
        //PRM_Name("eigen_cg", "Eigen's built-in cg solver"),
        //PRM_Name("pcg_direct_products", "Preconditioned CG - Normal"),
        PRM_Name(0)
    };
    static PRM_ChoiceList theSolverTypeMenu(PRM_CHOICELIST_SINGLE, theSolverTypeChoices);

    static PRM_Template theTemplates[] =
    {
        // Required Fluid Inputs
        PRM_Template(PRM_STRING,    1, &theVelocityName,            &theVelocityNameDefault),
        PRM_Template(PRM_STRING,    1, &theValidName,               &theValidNameDefault),
        PRM_Template(PRM_STRING,    1, &theViscosityName,           &theViscosityNameDefault),
        PRM_Template(PRM_STRING,    1, &theDensityName,             &theDensityNameDefault),
        PRM_Template(PRM_FLT,       1, &theMinDensityName,          PRMoneDefaults, 0, 0, 0, &PRM_SpareData::unitsDensity),
        PRM_Template(PRM_FLT,       1, &theMaxDensityName,          &theMaxDensityDefault, 0, 0, 0, &PRM_SpareData::unitsDensity),
        PRM_Template(PRM_STRING,    1, &thePressureName,            &thePressureNameDefault),
        PRM_Template(PRM_STRING,    1, &theSurfaceName,             &theSurfaceNameDefault),
        PRM_Template(PRM_STRING,    1, &theCollisionName,           &theCollisionNameDefault),
        PRM_Template(PRM_STRING,    1, &theCollisionVelocityName,   &theCollisionVelocityNameDefault),
        // Destination Outputs
        //PRM_Template(PRM_STRING, 1, &thePrintedGeoName, &thePrintedGeoDefault),
        // Optional Inputs
        PRM_Template(PRM_STRING,    1, &theExportDataPrefixName,    &theExportDataPrefixDefault),
        PRM_Template(PRM_ORD,       1, &theMatrixSchemeName, PRMoneDefaults, &theMatrixSchemeMenu),
        PRM_Template(PRM_ORD,       1, &theSolverTypeName, PRMoneDefaults, &theSolverTypeMenu),
        PRM_Template(PRM_TOGGLE,    1, &theDoSolveName, PRMoneDefaults),
        PRM_Template(PRM_TOGGLE,    1, &theKeepNonConvergedResultsName, PRMoneDefaults),
        PRM_Template(PRM_TOGGLE,    1, &theExportMatricesName, PRMzeroDefaults),
        PRM_Template(PRM_TOGGLE,    1, &theExportComponentMatricesName, PRMzeroDefaults),
        PRM_Template(PRM_TOGGLE,    1, &theExportStatsName, PRMzeroDefaults),
        PRM_Template(PRM_TOGGLE,    1, &theUseWarmStartName, PRMoneDefaults),
        PRM_Template(PRM_FLT,       1, &theToleranceName, &theToleranceDefault),
        PRM_Template(PRM_FLT,       1, &theMaxSolverIterations, &theMaxSolverIterationsDefault),
        PRM_Template(PRM_TOGGLE,    1, &theUseInputSurfaceWeightsName, PRMoneDefaults),
        PRM_Template(PRM_STRING,    1, &theSurfaceWeightsName, &theSurfaceWeightsNameDefault),
        PRM_Template(PRM_TOGGLE,    1, &theUseInputCollisionWeights, PRMoneDefaults),
        PRM_Template(PRM_STRING,    1, &theCollisionWeightsName, &theCollisionWeightsNameDefault),
        PRM_Template(PRM_INT,       1, &theActiveLiquidBoundaryLayerSizeName, &theActiveLiquidBoundaryLayerSizeDefault),
        PRM_Template(PRM_INT,       1, &theActiveSolidBoundaryLayerSizeName, &theActiveSolidBoundaryLayerSizeDefault),
        PRM_Template(PRM_TOGGLE,    1, &theDoReducedRegionsName, PRMoneDefaults),
        PRM_Template(PRM_TOGGLE,    1, &theDoTileName, PRMoneDefaults),
        PRM_Template(PRM_INT,       1, &theTileSizeName, &theTileSizeDefault),
        PRM_Template(PRM_INT,       1, &theTilePaddingName, &theTilePaddingDefault),
        PRM_Template()
    };

    static SIM_DopDescription    theDopDescription(
        true,                               // Should we make a DOP?
        "hdk_polystokes",                   // Internal name of the DOP.
        "HDK Polynomial Stokes Solver",     // Label of the DOP
        "$OS",                              // Default data name
        classname(),                        // The type of this DOP, usually the class.
        theTemplates);                      // Template list for generating the DOP
    setGasDescription(theDopDescription);

    return &theDopDescription;
}

bool
HDK_PolyStokes::solveGasSubclass(SIM_Engine& engine,
    SIM_Object* obj,
    SIM_Time time,
    SIM_Time timestep)
{
    std::cout << "Starting solveGasSubclass, compiled " << compileDateTime << std::endl;

    myObj = obj;

    ////////////////////////////////////////////
    // Load in required fields.
    ////////////////////////////////////////////
    SIM_VectorField* velocityField = getVectorField(obj, GAS_NAME_VELOCITY);
    SIM_VectorField* validField = getVectorField(obj, "valid");
    const SIM_ScalarField* viscosityField = getScalarField(obj, "viscosity");
    const SIM_ScalarField* pressureField = getScalarField(obj, "pressure");
    const SIM_ScalarField* densityField = getScalarField(obj, "density");

    const SIM_ScalarField* surfaceField = getConstScalarField(obj, GAS_NAME_SURFACE);
    const SIM_VectorField* surfaceWeights = getVectorField(obj, "surfaceweights");

    const SIM_ScalarField* collisionField = getConstScalarField(obj, GAS_NAME_COLLISION);
    const SIM_VectorField* collisionWeights = getVectorField(obj, "collisionweights");
    const SIM_VectorField* collisionVelocityField = getConstVectorField(obj, GAS_NAME_COLLISIONVELOCITY);

    ////////////////////////////////////////////
    // Check that required fields exist and are structured correctly
    ////////////////////////////////////////////
    if (!velocityField)
    {
        addError(obj, SIM_MESSAGE, "Velocity field is missing.", UT_ERROR_WARNING);
        return false;
    }
    if (!velocityField->isFaceSampled())
    {
        addError(obj, SIM_MESSAGE, "Velocity field must be a staggered grid.", UT_ERROR_ABORT);
        return false;
    }

    if (!validField)
    {
        addError(obj, SIM_MESSAGE, "Valid field is missing.", UT_ERROR_ABORT);
        return false;
    }
    else if (!validField->isAligned(velocityField))
    {
        addError(obj, SIM_MESSAGE, "Valid field must align with the velocity field.", UT_ERROR_ABORT);
        return false;
    }

    if (!surfaceField)
    {
        addError(obj, SIM_MESSAGE, "Surface field is missing.", UT_ERROR_ABORT);
        return false;
    }
    if (!collisionField)
    {
        addError(obj, SIM_MESSAGE, "Collision field is missing.", UT_ERROR_ABORT);
        return false;
    }
    if (!viscosityField)
    {
        addError(obj, SIM_MESSAGE, "Viscosity field is missing.", UT_ERROR_ABORT);
        return false;
    }
    if (!pressureField)
    {
        addError(obj, SIM_MESSAGE, "Pressure field is missing.", UT_ERROR_ABORT);
        return false;
    }
    if (!densityField)
    {
        addError(obj, SIM_MESSAGE, "Density field is missing.", UT_ERROR_ABORT);
        return false;
    }
    const SIM_RawField& liquidDensity = *densityField->getField();
    fpreal32 constantLiquidDensity = 0.;
    if (!liquidDensity.field()->isConstant(&constantLiquidDensity))
    {
        addError(obj, SIM_MESSAGE, "Variable density is not currently supported", UT_ERROR_WARNING);
        return false;
    }
    if (!surfaceWeights && getUseInputSurfaceWeights())
    {
        addError(obj, SIM_MESSAGE, "User requested to use input surface weights but that field is missing.", UT_ERROR_ABORT);
        return false;
    }
    if (!collisionWeights && getUseInputCollisionWeights())
    {
        addError(obj, SIM_MESSAGE, "User requested to use input collision weights but that field is missing.", UT_ERROR_ABORT);
        return false;
    }

    ////////////////////////////////////////////
    // Get field configuration and other simulation parameters.
    ////////////////////////////////////////////
    const fpreal        dt = timestep;
    const fpreal        dx = velocityField->getVoxelSize(0).maxComponent();

    ////////////////////////////////////////////
    // Copy field data into RawFields
    ////////////////////////////////////////////
    const SIM_RawField& surfaceFieldData      = *surfaceField->getField();
    const SIM_RawField& viscosityFieldData    = *viscosityField->getField();
    const SIM_RawField& densityFieldData      = *densityField->getField();
    assert(&viscosityFieldData && &densityFieldData);

    ////////////////////////////////////////////
    // Setup the solver helper class
    ////////////////////////////////////////////
    Solver mySolver(*this
        , dx, dt
        , velocityField
        , collisionVelocityField
        , surfaceField
        , surfaceWeights
        , collisionField
        , collisionWeights
        , densityField
        , constantLiquidDensity
        , viscosityField);
    mySolver.setupClockStart();

    ////////////////////////////////////////////
    // Compute volume fraction weights
    ////////////////////////////////////////////
    {
        UT_PerfMonAutoSolveEvent event(this, "Build integration weights");
        //mySolver.buildIntegrationWeights();
        mySolver.buildIntegrationWeightsAlt();
    }

    ////////////////////////////////////////////
    // Classify cells
    ////////////////////////////////////////////
    {
        // sets labels to GENERICFLUID, SOLID, or UNSOLVED
        UT_PerfMonAutoSolveEvent event(this, "Classify cells");
        mySolver.classifyCells();
    }
    {
        if (mySolver.doReducedRegions())
        {
            // sets GENERICFLUID to ACTIVEFLUID or REDUCED
            UT_PerfMonAutoSolveEvent event(this, "Construct reduced regions");
            mySolver.constructReducedRegions();
        }
        else
        {
            // sets GENERICFLUID to ACTIVEFLUID
            UT_PerfMonAutoSolveEvent event(this, "Construct only active regions");
            mySolver.constructOnlyActiveRegions();
        }
    }
    {
        // sets labels to GENERICFLUID, SOLID, or UNSOLVED
        UT_PerfMonAutoSolveEvent event(this, "Classify faces");
        mySolver.classifyFaces();
    }
    {
        // sets labels to GENERICFLUID, SOLID, or UNSOLVED
        UT_PerfMonAutoSolveEvent event(this, "Classify edges");
        mySolver.classifyEdges();
    }

    ////////////////////////////////////////////
    // Build indices
    ////////////////////////////////////////////
    {
        if (mySolver.doReducedRegions())
        {
            UT_PerfMonAutoSolveEvent event(this, "Build reduced indices");
            mySolver.constructCenterReducedIndices();
            mySolver.constructFacesReducedIndices();
            mySolver.constructEdgesReducedIndices();
        }
    }
    {
        UT_PerfMonAutoSolveEvent event(this, "Build active indices");
        mySolver.constructCenterActiveIndices();
        mySolver.constructFacesActiveIndices();
        mySolver.constructEdgesActiveIndices();
    }

    ////////////////////////////////////////////
    // Build reduced region center of masses,
    //  least squares fits, and mass matrices
    ////////////////////////////////////////////
    {
        if (mySolver.doReducedRegions())
        {
            {
                UT_PerfMonAutoSolveEvent event(this, "Compute Center of Masses");
                mySolver.computeCenterOfMasses();
            }
            {
                UT_PerfMonAutoSolveEvent event(this, "Compute Least Squares Fits");
                mySolver.computeLeastSquaresFits();
            }
        }
    }

    ////////////////////////////////////////////
    // Build reduced region matrices
    // (mass matrices and viscous JD'DJ')
    ////////////////////////////////////////////
    {
        if (mySolver.doReducedRegions())
        {
            {
                UT_PerfMonAutoSolveEvent event(this, "Compute Reduced Region Mass Matrices");
                mySolver.computeReducedMassMatrices();
            }
            {
                // todo we can simplify this now
                // todo make the solver itself figure out which version it needs
                UT_PerfMonAutoSolveEvent event(this, "Compute Reduced Region Viscosity Matrices (2JD'uDJ')");
                if (getMatrixScheme() == HDK_PolyStokes_Options::MatrixScheme::ALL_DOFS)
                    mySolver.computeReducedViscosityMatricesInteriorOnly();
                else if (getMatrixScheme() == HDK_PolyStokes_Options::MatrixScheme::PRESSURE_STRESS)
                    mySolver.computeReducedViscosityMatricesInteriorOnly();
                else if (getMatrixScheme() == HDK_PolyStokes_Options::MatrixScheme::PRESSURE_VELOCITY)
                    mySolver.computeReducedViscosityMatricesInteriorOnly();
                else
                    mySolver.computeReducedViscosityMatrices();
            }
        }
    }

    ////////////////////////////////////////////
    // Build matrix blocks
    ////////////////////////////////////////////
    {
        UT_PerfMonAutoSolveEvent event(this, "Construct Matrix Blocks");
        mySolver.constructMatrixBlocks();
    }

    ////////////////////////////////////////////
    // Setup guess vectors
    ////////////////////////////////////////////
    {
        UT_PerfMonAutoSolveEvent event(this, "DEBUG: Read In Known Solution for Warm Starting");
        mySolver.initializeGuessVectors();
        if (getUseWarmStart())
            mySolver.constructGuessVectors();
    }
    ////////////////////////////////////////////
    // Assemble system
    ////////////////////////////////////////////
    {
        UT_PerfMonAutoSolveEvent event(this, "Assemble System");
        mySolver.assemble();
    }
    mySolver.setupClockEnd();

    ////////////////////////////////////////////
    // Print data for display and debugging
    ////////////////////////////////////////////
    {
        UT_PerfMonAutoSolveEvent event(this, "Printing Grid Data");
        mySolver.printAllData();
    }

    {
        //if (getDoUsePreconditioner())
        {
            UT_PerfMonAutoSolveEvent event(this, "Construct Preconditioner");
            mySolver.constructPreconditioner();
        }
    }
    {
        if (getExportMatrices())
        {
            UT_PerfMonAutoSolveEvent event(this, "Exporting Matrix Data");
            mySolver.exportMatrices((std::string)getExportDataPrefix());
        }
    }
    {
        if (getExportComponentMatrices())
        {
            UT_PerfMonAutoSolveEvent event(this, "Exporting Component Matrix Data");
            mySolver.exportComponentMatrices((std::string)getExportDataPrefix());
        }
    }

    ////////////////////////////////////////////
    // Solve
    ////////////////////////////////////////////
    Solver::SolverResult solverResult = Solver::SolverResult::INCOMPLETE;
    {
        if (getDoSolve())
        {
            std::cout << "Starting solve..." << std::endl;
            
            UT_PerfMonAutoSolveEvent event(this, "Solve");
            solverResult = mySolver.solve();

            std::cout << "Done solve." << std::endl;
            mySolver.printStats();

        }
    }

    ////////////////////////////////////////////
    // Address Solver Errors
    ////////////////////////////////////////////
    {
        if (solverResult == Solver::SolverResult::UNSUPPORTED_SOLVER)
        {
            addError(obj, SIM_MESSAGE, "Unsupported Solver.", UT_ERROR_ABORT);
            return false;
        }
    }

    {
        //mySolver.extractResiduals();
        //mySolver.printResidualData();
    }

    {
        if (getExportMatrices())
        {
            UT_PerfMonAutoSolveEvent event(this, "Exporting Solution Data");
            mySolver.exportMatricesPostSolve((std::string)getExportDataPrefix());
        }
    }

    {
        if (getExportStats())
        {
            UT_PerfMonAutoSolveEvent event(this, "Exporting Solve Stats");
            mySolver.exportStats((std::string)getExportDataPrefix());
        }
    }

    ////////////////////////////////////////////
    // Push results out
    ////////////////////////////////////////////
    {
        UT_PerfMonAutoSolveEvent event(this, "Build Valid Field");
        mySolver.buildValidFaces(*validField);
    }
    {
        if (solverResult == Solver::SolverResult::SUCCESS || getKeepNonConvergedResults())
        {
            UT_PerfMonAutoSolveEvent event(this, "Update Reduced Velocities");

            // need to compute out velocities if we have a matrix scheme that eliminated those dofs
            if (getMatrixScheme() == HDK_PolyStokes_Options::MatrixScheme::PRESSURE_STRESS)
            {
                mySolver.recoverVelocityFromPressureStress();
            }

            for (int axis : {0, 1, 2})
            {
                mySolver.applySolutionToVelocity(
                    *velocityField->getField(axis),
                    *validField->getField(axis),
                    axis
                );
            }
        }
    }

    ////////////////////////////////////////////
    // Let Houdini know these fields were changed
    ////////////////////////////////////////////
    if (getDoSolve())
    {
        if (solverResult == Solver::SolverResult::SUCCESS || getKeepNonConvergedResults())
        {
            velocityField->pubHandleModification();
            validField->pubHandleModification();
        }
        else if (solverResult == Solver::SolverResult::NOCONVERGE && !getKeepNonConvergedResults())
        {
            addError(obj, SIM_MESSAGE, "Solver did not converge, exiting...", UT_ERROR_ABORT);
        }
        else
        {
            addError(obj, SIM_MESSAGE, "Solver failed, exiting...", UT_ERROR_ABORT);
        }
    }

    std::cout << std::endl;
    return solverResult == Solver::SolverResult::SUCCESS;
}
