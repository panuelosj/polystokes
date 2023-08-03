#pragma once

#include "units.h"

template<typename SIM_Field_T>
class WeightsSet
{
public:
	explicit WeightsSet(
		SIM_RawField _center,
		SIM_RawField _faceX,
		SIM_RawField _faceY,
		SIM_RawField _faceZ,
		SIM_RawField _edgeYZ,
		SIM_RawField _edgeXZ,
		SIM_RawField _edgeXY
	) :
		center_(_center)
		, faceX_(_faceX)
		, faceY_(_faceY)
		, faceZ_(_faceZ)
		, edgeYZ_(_edgeYZ)
		, edgeXZ_(_edgeXZ)
		, edgeXY_(_edgeXY)
	{};
	virtual ~WeightsSet();

	SIM_Field_T* center()
	{
		return &center_;
	}
	SIM_Field_T* face(int axis)
	{
		switch (axis)
		{
		case 0:
			return &faceX_;
		case 1:
			return &faceY_;
		case 2:
			return &faceZ_;
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}
	SIM_Field_T* edge(int axis)
	{
		switch (axis)
		{
		case 0:
			return &edgeYZ_;
		case 1:
			return &edgeXZ_;
		case 2:
			return &edgeXY_;
		default:
			assert(axis >= 0 && axis < 3);
			return NULL;
		}
	}

private:
	SIM_Field_T center_;
	SIM_Field_T faceX_, faceY_, faceZ_;
	SIM_Field_T edgeYZ_, edgeXZ_, edgeXZ_;
};

class GridClassifier
{
public:
	// Constructor and Destructor
	explicit GridClassifier();
	virtual ~GridClassifier();

private:
	void classifyCells();
	void constructReducedRegions();

	WeightsSet<SIM_RawField> liquidWeights, fluidWeights;
	WeightsSet<SIM_RawIndexField> labels;
	SIM_RawField mySurfaceFieldData;
	SIM_RawField myCollisionFieldData;
};