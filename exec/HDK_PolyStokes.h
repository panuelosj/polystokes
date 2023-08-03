/*
 * 
 */
#ifndef __HDK_PolyStokes_h__
#define __HDK_PolyStokes_h__

//#define CUDA_NO_HALF // cuda conflicts with houdini's half def

#include "units.h"
#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>
#include <SIM/SIM_FieldUtils.h>
#include <SIM/SIM_Object.h>

class HDK_PolyStokes : public GAS_SubSolver
{
public:
	// my solver
	class Solver;
	/// These macros are used to create the accessors
	/// getFieldDstName and getFieldSrcName functions we'll use
	/// to access our data options.
	GET_DATA_FUNC_E("matrixSetup",						MatrixScheme, HDK_PolyStokes_Options::MatrixScheme);
	GET_DATA_FUNC_E("solverType",						SolverType, HDK_PolyStokes_Options::SolverType);
	GET_DATA_FUNC_B("useInputSurfaceWeights",			UseInputSurfaceWeights);
	GET_DATA_FUNC_B("useInputCollisionWeights",			UseInputCollisionWeights);
	GET_DATA_FUNC_F("minDensity",						MinDensity);
	GET_DATA_FUNC_F("maxDensity",						MaxDensity);
	GET_DATA_FUNC_I("activeLiquidBoundaryLayerSize",	ActiveLiquidBoundaryLayerSize);
	GET_DATA_FUNC_I("activeSolidBoundaryLayerSize",		ActiveSolidBoundaryLayerSize);
	GET_DATA_FUNC_B("doReducedRegions",					DoReducedRegions);
	GET_DATA_FUNC_B("doTile",							DoTile);
	GET_DATA_FUNC_I("tileSize",							TileSize);
	GET_DATA_FUNC_I("tilePadding",						TilePadding); 
	GET_DATA_FUNC_F(SIM_NAME_TOLERANCE,					SolverTolerance);
	GET_DATA_FUNC_I("maxSolverIterations",				SolverMaxIterations);
	GET_DATA_FUNC_B("useWarmStart",						UseWarmStart);
	GET_DATA_FUNC_B("exportMatrices",					ExportMatrices);
	GET_DATA_FUNC_B("exportComponentMatrices",			ExportComponentMatrices);
	GET_DATA_FUNC_B("exportStats",						ExportStats);
	GET_DATA_FUNC_B("doSolve",							DoSolve);
	GET_DATA_FUNC_B("keepNonConvergedResults",			KeepNonConvergedResults);
	GET_DATA_FUNC_S("exportDataPrefix",					ExportDataPrefix);

protected:
	// Constructor and Destructor
	explicit	HDK_PolyStokes(const SIM_DataFactory* factory);
	virtual	   ~HDK_PolyStokes();

	/// The overloaded callback that GAS_SubSolver will invoke to
	/// perform our actual computation.  We are giving a single object
	/// at a time to work on.
	virtual bool solveGasSubclass(
		SIM_Engine& engine,
		SIM_Object* obj,
		SIM_Time time,
		SIM_Time timestep);

private:
	/// These macros are used to create the accessor functions
	/// to be able to set these data options.
	SET_DATA_FUNC_S("centerLabels", CenterLabels);
	SET_DATA_FUNC_S("faceXLabels", FaceXLabelsName);
	SET_DATA_FUNC_S("faceYLabels", FaceYLabelsName);
	SET_DATA_FUNC_S("faceZLabels", FaceZLabelsName);
	SET_DATA_FUNC_S("edgeXYLabels", EdgeXYLabelsName);
	SET_DATA_FUNC_S("edgeXZLabels", EdgeXZLabelsName);
	SET_DATA_FUNC_S("edgeYZLabels", EdgeYZLabelsName);

	SET_DATA_FUNC_S("centerReducedIndices", CenterReducedIndices);
	SET_DATA_FUNC_S("faceXReducedIndices", FaceXReducedIndicesName);
	SET_DATA_FUNC_S("faceYReducedIndices", FaceYReducedIndicesName);
	SET_DATA_FUNC_S("faceZReducedIndices", FaceZReducedIndicesName);
	SET_DATA_FUNC_S("edgeXYReducedIndices", EdgeXYReducedIndicesName);
	SET_DATA_FUNC_S("edgeXZReducedIndices", EdgeXZReducedIndicesName);
	SET_DATA_FUNC_S("edgeYZReducedIndices", EdgeYZReducedIndicesName);

	SET_DATA_FUNC_S("centerActiveIndices",CenterActiveIndices);
	SET_DATA_FUNC_S("faceXActiveIndices", FaceXActiveIndicesName);
	SET_DATA_FUNC_S("faceYActiveIndices", FaceYActiveIndicesName);
	SET_DATA_FUNC_S("faceZActiveIndices", FaceZActiveIndicesName);
	SET_DATA_FUNC_S("edgeXYActiveIndices", EdgeXYActiveIndicesName);
	SET_DATA_FUNC_S("edgeXZActiveIndices", EdgeXZActiveIndicesName);
	SET_DATA_FUNC_S("edgeYZActiveIndices", EdgeYZActiveIndicesName);

	SET_DATA_FUNC_S("centerLiquidWeights",	CenterLiquidWeightsName);
	SET_DATA_FUNC_S("faceXLiquidWeights",	FaceXLiquidWeightsName);
	SET_DATA_FUNC_S("faceYLiquidWeights",	FaceYLiquidWeightsName);
	SET_DATA_FUNC_S("faceZLiquidWeights",	FaceZLiquidWeightsName);
	SET_DATA_FUNC_S("edgeXYLiquidWeights",	EdgeXYLiquidWeightsName);
	SET_DATA_FUNC_S("edgeXZLiquidWeights",	EdgeXZLiquidWeightsName);
	SET_DATA_FUNC_S("edgeYZLiquidWeights",	EdgeYZLiquidWeightsName);

	SET_DATA_FUNC_S("centerFluidWeights",	CenterFluidWeightsName);
	SET_DATA_FUNC_S("faceXFluidWeights",	FaceXFluidWeightsName);
	SET_DATA_FUNC_S("faceYFluidWeights",	FaceYFluidWeightsName);
	SET_DATA_FUNC_S("faceZFluidWeights",	FaceZFluidWeightsName);
	SET_DATA_FUNC_S("edgeXYFluidWeights",	EdgeXYFluidWeightsName);
	SET_DATA_FUNC_S("edgeXZFluidWeights",	EdgeXZFluidWeightsName);
	SET_DATA_FUNC_S("edgeYZFluidWeights",	EdgeYZFluidWeightsName);

	SET_DATA_FUNC_S("reducedRegionCenterOfMass", ReducedRegionCenterOfMassName);

	// fields
	SIM_Object* myObj;

	/// We define this to be a DOP_Auto node which means we do not
	/// need to implement a DOP_Node derivative for this data.  Instead,
	/// this description is used to define the interface.
	static const SIM_DopDescription* getDopDescription();
	/// These macros are necessary to bind our node to the factory and
	/// ensure useful constants like BaseClass are defined.
	DECLARE_STANDARD_GETCASTTOTYPE();
	DECLARE_DATAFACTORY(HDK_PolyStokes,
		GAS_SubSolver,
		"HDK PolyStokes Solver",
		getDopDescription());

};
#endif