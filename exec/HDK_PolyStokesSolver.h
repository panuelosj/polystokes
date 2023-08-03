#ifndef __HDK_PolyStokes_Solver_h__
#define __HDK_PolyStokes_Solver_h__

#include "units.h"
#include "util.h"
//#define CUDA_NO_HALF // cuda conflicts with houdini's half def

#include "HDK_PolyStokes.h"
//#include <Eigen/IterativeLinearSolvers>
// solvers
#include <unsupported/Eigen/IterativeSolvers>
// preconditioner
#include "Preconditioner.h"
// applies the action of matrix A
#include "ApplyPressureStressMatrix.h"

// gas solver tools
#include <SIM/SIM_Object.h>
// multithreading
#include <UT/UT_ParallelUtil.h>
#include <tbb/tbb.h>

// exporting data 
#include <string>
#include <chrono>

class HDK_PolyStokes::Solver
{
public:
	// Constructor and Destructor
	// todo support variable density?
	explicit Solver(
		HDK_PolyStokes& _parent,
		const fpreal dx, const fpreal dt,
		SIM_VectorField* _velocityField,
		const SIM_VectorField* _collisionVelocityField,
		const SIM_ScalarField* _surfaceField,
		const SIM_VectorField* _surfaceWeights,
		const SIM_ScalarField* _collisionField,
		const SIM_VectorField* _collisionWeights,
		const SIM_ScalarField* _densityField,
		const fpreal _constantDensity,
		const SIM_ScalarField* _viscosityField);
	virtual ~Solver();

	enum class ValidFlag
	{
		VALID_FACE = 1,
		INVALID_FACE = 0
	};
	enum class SamplingType : int
	{
		CENTER = SIM_SAMPLE_CENTER,
		EDGEXY = SIM_SAMPLE_EDGEXY,
		EDGEXZ = SIM_SAMPLE_EDGEXZ,
		EDGEYZ = SIM_SAMPLE_EDGEYZ,
		FACEX = SIM_SAMPLE_FACEX,
		FACEY = SIM_SAMPLE_FACEY,
		FACEZ = SIM_SAMPLE_FACEZ
	};
	enum class SolverResult
	{
		UNSUPPORTED_SOLVER = -4,
		INCOMPLETE = -3,
		INVALID = -2,
		FAILED = -1,
		NOCONVERGE = 0,
		SUCCESS = 1,
		NOCHANGE = 2
	};
	enum MaterialLabels
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
	enum StressType
	{
		XX,
		YY,
		ZZ,
		YZ,
		XZ,
		XY
	};
	bool isStressType(StressType type) {
		return (type >= StressType::XX && type <= StressType::XY);
	}

	void buildIntegrationWeights();
	void buildIntegrationWeightsAlt();
	void computeIntegrationWeights(SIM_RawField& integrationWeights,
		const SIM_RawField& liquidSurface,
		const SIM_FieldSample sample,
		const int numberOfSamples);
	void computeSolidIntegrationWeights(SIM_RawField& solidIntegrationWeights,
		const SIM_RawField& liquidSurface,
		const SIM_RawField& solidSurface,
		const SIM_FieldSample sample,
		const int numberOfSamples,
		const fpreal extrapolation);
	void classifyCells();
	void classifyCellsAlt();
	void constructReducedRegions();
	void constructOnlyActiveRegions();
	void classifyFaces();
	void classifyEdges();

	void constructCenterReducedIndices();
	void constructFacesReducedIndices();
	void constructEdgesReducedIndices();

	void constructCenterActiveIndices();
	void constructFacesActiveIndices();
	void constructEdgesActiveIndices();

	void computeCenterOfMasses();
	void computeLeastSquaresFits();
	void computeReducedMassMatrices();
	void computeReducedViscosityMatrices();
	void computeReducedViscosityMatricesInteriorOnly();

	void constructMatrixBlocks();

	void assemble();
	void assembleSystem();
	void assembleSystemAlt();
	void assembleSystemExplicitInternalStresses();
	void assembleSystemVelocityPressure();
	void assembleSystemPressureStress();
	void assembleSystemPressureStressFactored();
	void exportMatrices(std::string prefix);
	void exportComponentMatrices(std::string prefix);
	void exportMatricesPostSolve(std::string prefix);
	void exportStats(std::string prefix);
	void printStats();
	void initializeGuessVectors();
	void constructGuessVectors();
	void readInWarmStart();

	// preconditioners if we need to use solvePCG()
	void constructPreconditioner();
	void constructPreconditionerGSsmoother();
	void constructPreconditionerIdentity();
	void constructPreconditionerDiagonal();
	void constructPreconditionerEq6();
	void constructPreconditionerEq14();

	// solvers
	SolverResult solve();
	bool checkSolverMatrixCompatibility();
	// general solvers
	SolverResult solvePCG();
	SolverResult solveEigenCG();
	// specific solvers
	SolverResult solveSPDwithMatrixVectorPCG();
	SolverResult solveGSSmoother();
	SolverResult solveGSSmootherNew();
	void stepGSSmootherUniform(Vector& soln_us, Vector& soln_vs, Vector& soln_ps);
	void stepGSSmootherReduced(Vector& soln_us, Vector& soln_vs, Vector& soln_ps);
	
	void extractResiduals();
	void printResidualData();
	void writeVectorToField(Vector data, SIM_RawField& field, SamplingType type);

	void buildValidFaces(SIM_VectorField& validFaces);
	void recoverVelocityFromPressureStress();
	void applySolutionToVelocity(SIM_RawField& velocity, SIM_RawField& validFaces, const int axis);

	// output stuff
	void printAllData();
	void printData(SIM_RawField& data, SamplingType sampling, const char* name);
	void printIndexData(SIM_RawIndexField& data, SamplingType sampling, const char* name);
	void printArrayData(UT_Array<UT_Vector3T<SolveReal>>& data, const char* name);

	// timing
	void setupClockStart();
	void setupClockEnd();

	// variable passthrough
	bool doReducedRegions() {
		return myDoReducedRegions;
		//return false;
	}

protected:
	UT_Vector3 SamplingOffset(SamplingType type)
	{
		switch (type)
		{
		case SamplingType::CENTER:
			return UT_Vector3(0.5, 0.5, 0.5);
			break;
		case SamplingType::EDGEXY:
			return UT_Vector3(0., 0., 0.5);
			break;
		case SamplingType::EDGEXZ:
			return UT_Vector3(0., 0.5, 0.);
			break;
		case SamplingType::EDGEYZ:
			return UT_Vector3(0.5, 0., 0.);
			break;
		case SamplingType::FACEX:
			return UT_Vector3(0., 0.5, 0.5);
			break;
		case SamplingType::FACEY:
			return UT_Vector3(0.5, 0., 0.5);
			break;
		case SamplingType::FACEZ:
			return UT_Vector3(0.5, 0.5, 0.);
			break;
		default:
			return UT_Vector3(0., 0., 0.);
			break;
		}
	}

private:
	// integration constants
	static constexpr fpreal MINWEIGHT = 0.1;
	static constexpr int NSAMPLES = 2;

	// inputs
	HDK_PolyStokes& myParent;
	SIM_VectorField* myVelocityField;
	const SIM_VectorField* myCollisionVelocityField;
	const bool myUseSurfaceWeights;
	const SIM_ScalarField* mySurfaceField;
	const SIM_VectorField* mySurfaceWeights;
	const bool myUseCollisionWeights;
	const SIM_ScalarField* myCollisionField;
	const SIM_VectorField* myCollisionWeights;

	const SIM_RawField* mySurfaceFieldData;
	const SIM_RawField* myCollisionFieldData;
	const SIM_RawField* myViscosityFieldData;
	const SIM_RawField* myDensityFieldData;
	const fpreal myConstantDensity;
	const fpreal myMinDensity;
	const fpreal myMaxDensity;

	// solver options
	const HDK_PolyStokes_Options::MatrixScheme myMatrixScheme;
	const HDK_PolyStokes_Options::SolverType mySolverType;
	const fpreal mySolverTolerance;
	const int mySolverMaxIterations;
	fpreal setupCPUTime = -1., setupWallclockTime = -1.;
	fpreal solveCPUTime = -1., solveWallclockTime = -1., solveError = -1.;
	int solveIterations = -1;

	const fpreal dx, invDx, invDx2;
	const fpreal dt, invDt;
	const UT_Vector3 size;
	const UT_Vector3 orig;
	const UT_Vector3 resolution;
	const exint nx;
	const exint ny;
	const exint nz;

	SolverResult mySolverResult = SolverResult::INCOMPLETE;

	// clock
	std::clock_t setupCPUClockStart, setupCPUClockEnd;
	std::chrono::time_point<std::chrono::high_resolution_clock> setupWallClockStart, setupWallClockEnd;

	// number of active grid points
	exint nCenter = 0;
	exint nFaceX = 0, nFaceY = 0, nFaceZ = 0;
	exint nEdgeXY = 0, nEdgeXZ = 0, nEdgeYZ = 0;
	// number of dofs
	exint nActiveVs = 0;
	exint nActiveVx = 0, nActiveVy = 0, nActiveVz = 0;
	exint nReducedVs = 0;
	exint nPressures = 0;
	exint nStresses = 0;
	exint nReducedStresses = 0;
	exint nTxx = 0, nTyy = 0, nTzz = 0;
	exint nTxy = 0, nTxz = 0, nTyz = 0;
	exint nTotalDOFs = 0, nSystemSize = 0;

	int myThreadCount, myGrainSize;
	int myLiquidBoundaryLayerSize, mySolidBoundaryLayerSize;
	bool myDoReducedRegions, myDoTile;
	int myTileSize, myTilePadding;
	exint myInteriorRegionCount = 0;

	// out of bounds checks
	bool c_oob(int i, int j, int k) const {
		return i < 0 || i > nx - 1 || j < 0 || j > ny - 1 || k < 0 || k > nz - 1;
	}
	bool tyz_oob(int i, int j, int k) const {
		return i < 0 || i > nx - 1 || j < 0 || j > ny || k < 0 || k > nz;
	}
	bool txz_oob(int i, int j, int k) const {
		return i < 0 || i > nx || j < 0 || j > ny - 1 || k < 0 || k > nz;
	}
	bool txy_oob(int i, int j, int k) const {
		return i < 0 || i > nx || j < 0 || j > ny || k < 0 || k > nz - 1;
	}
	bool x_oob(int i, int j, int k) const {
		return i < 0 || i > nx || j < 0 || j > ny - 1 || k < 0 || k > nz - 1;
	}
	bool y_oob(int i, int j, int k) const {
		return i < 0 || i > nx - 1 || j < 0 || j > ny || k < 0 || k > nz - 1;
	}
	bool z_oob(int i, int j, int k) const {
		return i < 0 || i > nx - 1 || j < 0 || j > ny - 1 || k < 0 || k > nz;
	}

	// computed in buildIntegrationWeights()
	std::vector<SIM_RawField*> liquidWeights;
	std::vector<SIM_RawField*> fluidWeights;
	SIM_RawField centerLiquidWeights, faceXLiquidWeights, faceYLiquidWeights, faceZLiquidWeights;
	SIM_RawField edgeXYLiquidWeights, edgeXZLiquidWeights, edgeYZLiquidWeights;
	SIM_RawField centerFluidWeights, faceXFluidWeights, faceYFluidWeights, faceZFluidWeights;
	SIM_RawField edgeXYFluidWeights, edgeXZFluidWeights, edgeYZFluidWeights;
	
	// computed in extractResiduals()
	SIM_RawField pressureResidual, velXResidual, velYResidual, velZResidual;

	// material labels
	SIM_RawIndexField centerLabels, faceXLabels, faceYLabels, faceZLabels;
	SIM_RawIndexField edgeXYLabels, edgeXZLabels, edgeYZLabels;
	// uniform indices
	SIM_RawIndexField centerActiveIndices, faceXActiveIndices, faceYActiveIndices, faceZActiveIndices;
	SIM_RawIndexField edgeXYActiveIndices, edgeXZActiveIndices, edgeYZActiveIndices;
	// reduced indices
	SIM_RawIndexField centerReducedIndices, faceXReducedIndices, faceYReducedIndices, faceZReducedIndices;
	SIM_RawIndexField edgeXYReducedIndices, edgeXZReducedIndices, edgeYZReducedIndices;

	// matrix blocks
	// diagonal blocks
	UT_Array<ReducedMatrix> reducedMassMatrices, reducedViscosityMatrices;
	SparseMatrix Mc_Matrix, McInv_Matrix;	// nActiveVs x nActiveVs
	SparseMatrix Mr_Matrix, MrInv_Matrix;	// nReducedVs x nReducedVs
	SparseMatrix JDtuDJ_Matrix;				// nReducedVs x nReducedVs
	SparseMatrix Mr_plus_2JDtuDJ_Matrix;	// nReducedVs x nReducedVs
	SparseMatrix Inv_Mr_plus_2JDtuDJ_Matrix;// nReducedVs x nReducedVs
	SparseMatrix B_Matrix;					// nReducedVs x nReducedVs	// todo remove this after debugging
	SparseMatrix BInv_Matrix;				// nReducedVs x nReducedVs
	SparseMatrix uInv_Matrix, u_Matrix;		// nStresses x nStresses
	SparseMatrix uInvRed_Matrix, uRed_Matrix;	// nReducedStresses x nReducedStresses
	SparseMatrix V_Matrix, VJt_Matrix, JVJt_Matrix;
	// off diagonal blocks
	SparseMatrix G_Matrix;					// nActiveVs x nPressures
	SparseMatrix Dt_Matrix;					// nActiveVs x nStresses
	SparseMatrix JG_Matrix;					// nReducedVs x nPressures
	SparseMatrix JDt_Matrix;				// nReducedVs x nStresses
	SparseMatrix JDtRed_Matrix;				// nReducedVs x nReducedStresses
	// rhs vectors
	Vector activeRHSVector;					// nActiveVs
	Vector reducedRHSVector;				// nReducedVs
	Vector pressureRHSVector;				// nPressures
	Vector stressRHSVector;					// nStresses
	SparseMatrix activeRHSVectorSparse, pressureRHSVectorSparse, stressRHSVectorSparse;
	// old velocities
	SparseMatrix oldActiveVs;				// nActiveVs
	SparseMatrix oldReducedVs;				// nReducedVs
	// guess vectors
	Vector activeGuessVector;				// nActiveVs
	Vector reducedGuessVector;				// nReducedVs
	Vector pressureGuessVector;				// nPressures
	Vector stressGuessVector;				// nStresses
	// system
	Preconditioner* presolver;
	SparseMatrix A;
	Vector b, guessVector, solutionVector;
	Vector residual;
	ApplyPressureStressMatrix<>* applyMatrix;


	SamplingType faceAxisToSamplingType(uint axis)
	{
		return (SamplingType)(axis + 4);
	}
	uint samplingTypeToFaceAxis(SamplingType type)
	{
		switch (type)
		{
		case SamplingType::FACEX:
			return 0;
		case SamplingType::FACEY:
			return 1;
		case SamplingType::FACEZ:
			return 2;
		default:
			assert(type == SamplingType::FACEX || type == SamplingType::FACEY || type == SamplingType::FACEZ);
			return -1;
		}
	}
	uint samplingTypeToEdgeAxis(SamplingType type)
	{
		switch (type)
		{
		case SamplingType::EDGEXY:
			return 2;
		case SamplingType::EDGEXZ:
			return 1;
		case SamplingType::EDGEYZ:
			return 0;
		default:
			assert(type == SamplingType::EDGEXY || type == SamplingType::EDGEXZ || type == SamplingType::EDGEYZ);
			return -1;
		}
	}

	Index gridLocationToActiveIndex(UT_Vector3I coord, SamplingType type)
	{
		using SIM::FieldUtils::getFieldValue;

		switch (type)
		{
		case SamplingType::EDGEXY:
		case SamplingType::EDGEXZ:
		case SamplingType::EDGEYZ:
			return getFieldValue(*edgeActiveIndices(samplingTypeToEdgeAxis(type)), coord);
			break;
		case SamplingType::FACEX:
		case SamplingType::FACEY:
		case SamplingType::FACEZ:
			return getFieldValue(*faceActiveIndices(samplingTypeToFaceAxis(type)), coord);
			break;
		case SamplingType::CENTER:
		default:
			return getFieldValue(centerActiveIndices, coord);
			break;
		}
	}

	SIM_RawField* faceFluidWeights(int axis)
	{
		switch (axis)
		{
		case 0:
			return &faceXFluidWeights;
		case 1:
			return &faceYFluidWeights;
		case 2:
			return &faceZFluidWeights;
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	SIM_RawField* faceLiquidWeights(int axis)
	{
		switch (axis)
		{
		case 0:
			return &faceXLiquidWeights;
		case 1:
			return &faceYLiquidWeights;
		case 2:
			return &faceZLiquidWeights;
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	SIM_RawField* edgeFluidWeights(int axis)
	{
		switch (axis)
		{
		case 0:
			return &edgeYZFluidWeights;
		case 1:
			return &edgeXZFluidWeights;
		case 2:
			return &edgeXYFluidWeights;
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	SIM_RawField* edgeLiquidWeights(int axis)
	{
		switch (axis)
		{
		case 0:
			return &edgeYZLiquidWeights;
		case 1:
			return &edgeXZLiquidWeights;
		case 2:
			return &edgeXYLiquidWeights;
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	SIM_RawIndexField* faceLabels(int axis)
	{
		switch (axis)
		{
		case 0:
			return &faceXLabels;
		case 1:
			return &faceYLabels;
		case 2:
			return &faceZLabels;
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	SIM_RawIndexField* faceActiveIndices(int axis)
	{
		switch (axis)
		{
		case 0:
			return &faceXActiveIndices;
		case 1:
			return &faceYActiveIndices;
		case 2:
			return &faceZActiveIndices;
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	SIM_RawIndexField* faceReducedIndices(int axis)
	{
		switch (axis)
		{
		case 0:
			return &faceXReducedIndices;
		case 1:
			return &faceYReducedIndices;
		case 2:
			return &faceZReducedIndices;
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	SIM_RawIndexField* edgeLabels(int axis)
	{
		switch (axis)
		{
		case 0:
			return &edgeYZLabels;
		case 1:
			return &edgeXZLabels;
		case 2:
			return &edgeXYLabels;
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	SIM_RawIndexField* edgeActiveIndices(int axis)
	{
		switch (axis)
		{
		case 0:
			return &edgeYZActiveIndices;
		case 1:
			return &edgeXZActiveIndices;
		case 2:
			return &edgeXYActiveIndices;
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	SIM_RawIndexField* edgeReducedIndices(int axis)
	{
		switch (axis)
		{
		case 0:
			return &edgeYZReducedIndices;
		case 1:
			return &edgeXZReducedIndices;
		case 2:
			return &edgeXYReducedIndices;
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	exint stressDOF(exint index, StressType type)
	{
		switch (type)
		{
		case StressType::XX:
			return index + 0;
		case StressType::YY:
			return index + nTxx;
		case StressType::ZZ:
			return index + nTxx + nTyy;
		case StressType::YZ:
			return index + nTxx + nTyy + nTzz;
		case StressType::XZ:
			return index + nTxx + nTyy + nTzz + nTyz;
		case StressType::XY:
			return index + nTxx + nTyy + nTzz + nTyz + nTxz;
		default:
			assert(isStressType(type));
			return NULL;
		}
	}
	exint reducedStressDOF(exint index, StressType type)
	{
		switch (type)
		{
		case StressType::XX:
			return index + myInteriorRegionCount * 0.;
		case StressType::YY:
			return index + myInteriorRegionCount * 1.;
		case StressType::ZZ:
			return index + myInteriorRegionCount * 2.;
		case StressType::YZ:
			return index + myInteriorRegionCount * 3.;
		case StressType::XZ:
			return index + myInteriorRegionCount * 4.;
		case StressType::XY:
			return index + myInteriorRegionCount * 5.;
		default:
			assert(isStressType(type));
			return NULL;
		}
	}
	exint faceVelocityDOF(exint index, int axis)
	{
		switch (axis)
		{
		case 0:
			return index;
		case 1:
			return index + nFaceX;
		case 2:
			return index + nFaceX + nFaceY;
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	exint centerStressDOF(exint index, int axis)
	{
		switch (axis)
		{
		case 0:
			return stressDOF(index, StressType::XX);
		case 1:
			return stressDOF(index, StressType::YY);
		case 2:
			return stressDOF(index, StressType::ZZ);
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	exint reducedCenterStressDOF(exint index, int axis)
	{
		switch (axis)
		{
		case 0:
			return reducedStressDOF(index, StressType::XX);
		case 1:
			return reducedStressDOF(index, StressType::YY);
		case 2:
			return reducedStressDOF(index, StressType::ZZ);
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	exint edgeStressDOF(exint index, int axis)
	{
		switch (axis)
		{
		case 0:
			return stressDOF(index, StressType::YZ);
		case 1:
			return stressDOF(index, StressType::XZ);
		case 2:
			return stressDOF(index, StressType::XY);
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	exint reducedEdgeStressDOF(exint index, int axis)
	{
		switch (axis)
		{
		case 0:
			return reducedStressDOF(index, StressType::YZ);
		case 1:
			return reducedStressDOF(index, StressType::XZ);
		case 2:
			return reducedStressDOF(index, StressType::XY);
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}

	// interior region data
	UT_Array<UT_Vector3T<SolveReal>> reducedRegionCOM;
	UT_Array<ColumnVector> reducedRegionBestFitVectors;

	bool isActive(exint label) {
		return (MaterialLabels)label == MaterialLabels::ACTIVEFLUID || (MaterialLabels)label == MaterialLabels::BOUNDARY;
	}
	bool isReduced(exint label) {
		return (MaterialLabels)label == MaterialLabels::REDUCED || (MaterialLabels)label == MaterialLabels::BOUNDARY;
	}
	bool isReducedButNotBoundary(exint label) {
		return (MaterialLabels)label == MaterialLabels::REDUCED;
	}
	bool isSolved(exint label) {
		return (MaterialLabels)label == MaterialLabels::GENERICFLUID
			|| (MaterialLabels)label == MaterialLabels::ACTIVEFLUID
			|| (MaterialLabels)label == MaterialLabels::REDUCED
			|| (MaterialLabels)label == MaterialLabels::BOUNDARY;
	}

	// called in constructReducedRegions()
	void constructAirBoundaryLayer();
	void buildInitialAirBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>>& parallelActiveCenterList);
	void buildNextLiquidBoundaryLayer(
		UT_Array<UT_Array<UT_Vector3I> >& parallelNewActiveCellLayer,
		const UT_Array<UT_Vector3I>& oldActiveCellLayer);
	void constructSolidBoundaryLayer();
	void buildInitialSolidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>>& parallelActiveCenterList);
	void buildNextSolidBoundaryLayer(
		UT_Array<UT_Array<UT_Vector3I> >& parallelNewActiveCellLayer,
		const UT_Array<UT_Vector3I>& oldActiveCellLayer,
		const SIM_RawIndexField& visitedCells);
	void constructTiles();

	//
	THREADED_METHOD(HDK_PolyStokes, centerLabels.shouldMultiThread(),
		classifyCellsAxis)
		void classifyCellsAxisPartial(
			const UT_JobInfo& info);
	MaterialLabels findCenterLabel(uint i, uint j, uint k);

	// called in classifyFaces()
	THREADED_METHOD2(HDK_PolyStokes, faceLabels.shouldMultiThread(),
		classifyFaceAxis,
		SIM_RawIndexField&, faceLabels,
		SamplingType, type)
		void classifyFaceAxisPartial(
			SIM_RawIndexField& faceLabels,
			SamplingType type,
			const UT_JobInfo& info);
	MaterialLabels findFaceLabelFromCenter(uint i, uint j, uint k, SamplingType type);
	MaterialLabels findFaceLabelFromCenterAlt(uint i, uint j, uint k, SamplingType type);

	// called in classifyEdges()
	THREADED_METHOD2(HDK_PolyStokes, edgeLabels.shouldMultiThread(),
		classifyEdgeAxis,
		SIM_RawIndexField&, edgeLabels,
		SamplingType, type)
		void classifyEdgeAxisPartial(
			SIM_RawIndexField& edgeLabels,
			SamplingType type,
			const UT_JobInfo& info);
	MaterialLabels findEdgeLabelFromFace(uint i, uint j, uint k, SamplingType type);
	MaterialLabels findEdgeLabelFromFaceAlt(uint i, uint j, uint k, SamplingType type);

	// called in constructCenterReducedIndices()
	void fixReducedRegionBoundaries();
	void shrinkToBoxReducedRegions(
		UT_Array<UT_Vector3I>& interiorRegionBBMin,
		UT_Array<UT_Vector3I>& interiorRegionBBMax);
	void fixSmallReducedRegions();
	void buildInteriorBoundingBoxes(UT_Array<UT_Array<UT_Vector3I> >& parallelInteriorRegionBBMin,
		UT_Array<UT_Array<UT_Vector3I> >& parallelInteriorRegionBBMax);

	// called in constructFacesReducedIndices()
	THREADED_METHOD3(HDK_PolyStokes, faceReducedIndexField.shouldMultiThread(),
		constructFaceAxisReducedIndices,
		SIM_RawIndexField&, faceReducedIndexField,
		SIM_RawIndexField&, faceLabels,
		SamplingType, type)
		void constructFaceAxisReducedIndicesPartial(
			SIM_RawIndexField& faceReducedIndexField,
			SIM_RawIndexField& faceLabels,
			SamplingType type,
			const UT_JobInfo& info);
	exint findFaceReducedIndexFromCenter(uint i, uint j, uint k, SamplingType type);

	// called in constructEdgesReducedIndices()
	THREADED_METHOD3(HDK_PolyStokes, edgeReducedIndexField.shouldMultiThread(),
		constructEdgeAxisReducedIndices,
		SIM_RawIndexField&, edgeReducedIndexField,
		SIM_RawIndexField&, edgeLabels,
		SamplingType, type)
		void constructEdgeAxisReducedIndicesPartial(
			SIM_RawIndexField& edgeReducedIndexField,
			SIM_RawIndexField& edgeLabels,
			SamplingType type,
			const UT_JobInfo& info);
	exint findEdgeReducedIndexFromFaces(uint i, uint j, uint k, SamplingType type);

	// called in construct___ActiveIndices()
	exint serialAssignFieldIndices(
		SIM_RawIndexField& indexField,
		SIM_RawIndexField& labels,
		exint startIndex);

	// called in computeCenterOfMasses()
	void buildReducedRegionCOM(
		UT_Array<UT_Array<UT_Vector3T<SolveReal>>>& parallelInteriorRegionCOM,
		UT_Array<UT_Array<SolveReal>>& parallelInteriorRegionCellCount);

	
	// called in computeLeastSquaresFit()
	void buildInteriorBestFitSystems(
		UT_Array<UT_Array<ReducedMatrix>>& parallelBestFitMatrix,
		UT_Array<UT_Array<ColumnVector>>& parallelBestFitRHS);
		
	// called in computeReducedMassMatrices()
	void buildReducedMassMatrixSystems(UT_Array<UT_Array<ReducedMatrix>>& parallelReducedMassMatrix);

	// called in computeReducedViscosityMatrices()
	void buildReducedViscosityMatrixSystems(UT_Array<UT_Array<ReducedMatrix>>& parallelReducedViscosityMatrix);
	void buildReducedViscosityMatrixSystemsInteriorOnly(UT_Array<UT_Array<ReducedMatrix>>& parallelReducedViscosityMatrix);

	// called in constructMatrixBlocks()
	void buildMatrixBlocksByTriplets(
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelActiveMassMatrixElements,
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelActiveMassInverseMatrixElements,
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelActiveStressInverseMatrixElements,
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelActiveStressMatrixElements,
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelReducedStressInverseMatrixElements,
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelReducedStressMatrixElements,
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelActiveGradientMatrixElements,
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelReducedGradientMatrixElements,
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelActiveDivergenceMatrixElements,
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelReducedDivergenceMatrixElements,
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelReducedInternalDivergenceMatrixElements,
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelActiveRHSElements,
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelPressureRHSElements,
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelStressRHSElements,
		std::vector<std::vector<Eigen::Triplet<SolveReal>>>& parallelActiveOldVsElements
	);

	// called in assembleSystem()
	void assembleReducedMassBlock();
	void assembleReducedMassBlockInverse();
	void assembleReducedViscosityBlock();
	void assembleReducedCombinedBlock();
	void assembleReducedInvertedBlock();
	void assembleReducedInvertedBlockIncludingBoundaryStresses();
	void assembleReducedUninvertedBlockIncludingBoundaryStresses();
	void assembleReducedRHSVector();
	void assembleVMatrices();

	// reduced conversions
	ColumnVector buildConversionCoefficients(UT_Vector3T<SolveReal> offset, int axis);

	// utils
	fpreal getLocalDensity(UT_Vector3I face);
	fpreal getLocalViscosity(UT_Vector3 point);
	void copyMaterialLabel(
		SIM_RawIndexField& source,
		SIM_RawIndexField& dest,
		const exint searchValue);
	void overwriteIndices(
		SIM_RawIndexField& indices,
		const exint searchValue,
		const exint replaceValue);
	template<typename Grid> void initField(
		Grid& weights,
		const SIM_RawField& surface,
		const SIM_FieldSample sample);
	template<typename Grid> void uncompressTiles(
		Grid& grid, 
		const UT_Array<bool>& isTileOccupiedList);
	void setActiveLayerCells(const UT_Array<UT_Vector3I>& activeCellLayer);
	void findOccupiedIndexTiles(
		UT_Array<bool>& isTileOccupiedList,
		const UT_Array<UT_Vector3I>& indexCellList,
		const SIM_RawIndexField& indexCellLabels);
	void remapInteriorRegions(
		const UT_Array<bool>& doRemoveRegion,
		const UT_Array<exint>& remapRegion);
};

#endif