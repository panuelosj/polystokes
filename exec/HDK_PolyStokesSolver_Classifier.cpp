#include "HDK_PolyStokesSolver.h"
#include <SIM/SIM_VolumetricConnectedComponentBuilder.h>

void
HDK_PolyStokes::Solver::buildValidFaces(SIM_VectorField& validFaces)
{
    using namespace SIM::FieldUtils;

    UT_Array<bool> isTileOccupiedList;
    // Uncompress valid face tiles

    for (int axis : {0, 1, 2})
    {
        validFaces.getField(axis)->makeConstant((float)ValidFlag::VALID_FACE);
        isTileOccupiedList.clear();
        isTileOccupiedList.setSize(validFaces.getField(axis)->field()->numTiles());
        isTileOccupiedList.constant(true);

        /*	HDK_Utilities::findOccupiedFaceTiles(isTileOccupiedList,
                                *validFaces.getField(axis),
                                faceVolumes[axis],
                                axis);
                                */
        uncompressTiles(*validFaces.getField(axis), isTileOccupiedList);

        UTparallelForEachNumber(validFaces.getField(axis)->field()->numTiles(), [&](const UT_BlockedRange<int>& range)
            {
                UT_VoxelArrayIteratorF vit;
                vit.setConstArray(validFaces.getField(axis)->fieldNC());
                vit.setCompressOnExit(true);

                for (exint i = range.begin(); i != range.end(); ++i)
                {
                    vit.myTileStart = i;
                    vit.myTileEnd = i + 1;
                    vit.rewind();

                    // Tiles for the index grid must have already been uncompressed
                    if (!vit.isTileConstant())
                    {
                        for (; !vit.atEnd(); vit.advance())
                        {
                            UT_Vector3I face(vit.x(), vit.y(), vit.z());

                            if (getFieldValue(*faceLabels(axis), face) == MaterialLabels::UNSOLVED || getFieldValue(*faceLabels(axis), face) == MaterialLabels::UNASSIGNED)
                                vit.setValue((float)ValidFlag::INVALID_FACE);
                            else
                                vit.setValue((float)ValidFlag::VALID_FACE);
                        }
                    }
                }
            });
    }
}

void
HDK_PolyStokes::Solver::classifyCells()
{
    // in: centerLabels = MaterialLabels::UNASSIGNED
    // out: centerLabels = MaterialLabels::{SOLID, GENERICFLUID, UNSOLVED}
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;
    using SIM::FieldUtils::cellToFaceMap;

    UT_Interrupt* boss = UTgetInterrupt();

    UTparallelForEachNumber(mySurfaceFieldData->field()->numTiles(), [&](const UT_BlockedRange<int>& range)
        {
            UT_VoxelArrayIteratorF vit;
            vit.setConstArray(mySurfaceFieldData->field());
            UT_VoxelTileIteratorF vitt;

            if (boss->opInterrupt())
                return;

            for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
            {
                vit.myTileStart = tileNumber;
                vit.myTileEnd = tileNumber + 1;
                vit.rewind();

                if (!vit.atEnd())
                {
                    vitt.setTile(vit);
                    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                    {
                        bool isInSolve = false;
                        bool isInFluid = true;

                        UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

                        // any cell with nonzero weights are solved
                        if (getFieldValue(centerLiquidWeights, cell) > 0.)
                            isInSolve = true;

                        // if any adjacent faces are inSolve, then cell should also be inSolve
                        if (!isInSolve)
                            for (int axis : {0, 1, 2})
                                for (int direction : {0, 1})
                                {
                                    UT_Vector3I face = cellToFaceMap(cell, axis, direction);
                                    if (getFieldValue(*faceLiquidWeights(axis), face) > 0.)
                                        isInSolve = true;
                                }

                        if (getFieldValue(centerFluidWeights, cell) == 0.)
                            isInFluid = false;

                        if (isInSolve)
                        {
                            if (!isInFluid)
                                setFieldValue(centerLabels, cell, MaterialLabels::SOLID);
                            else
                                setFieldValue(centerLabels, cell, MaterialLabels::GENERICFLUID);
                        }
                        else
                        {
                            setFieldValue(centerLabels, cell, MaterialLabels::UNSOLVED);
                        }
                    }
                }
            }
        });

    // we just wrote into this field, so make sure to compress the data in case entire tiles are constant
    // // in progress this might be neededj
    //centerLabels.fieldNC()->collapseAllTiles();
}

void
HDK_PolyStokes::Solver::classifyCellsAlt()
{
    // unused
    classifyCellsAxis();
}

void
HDK_PolyStokes::Solver::classifyCellsAxisPartial(
    const UT_JobInfo& info)
{
    using SIM::FieldUtils::getFieldValue;

    UT_VoxelArrayIteratorI vit(centerLabels.fieldNC());
    vit.splitByTile(info);

    for (vit.rewind(); !vit.atEnd(); vit.advance())
    {
        int i = vit.x(), j = vit.y(), k = vit.z();
        UT_Vector3I edge(vit.x(), vit.y(), vit.z());

        MaterialLabels label = MaterialLabels::UNASSIGNED;
        label = findCenterLabel(i, j, k);
        centerLabels.fieldNC()->setValue(UT_Vector3I(i, j, k), label);
    }
}

HDK_PolyStokes::Solver::MaterialLabels
HDK_PolyStokes::Solver::findCenterLabel(uint i, uint j, uint k)
{
    bool insystem = false;
    insystem = (!c_oob(i, j, k) && centerFluidWeights.fieldNC()->getValue(i, j, k));

    if (!insystem)
    {
        return MaterialLabels::UNSOLVED;
    }

    insystem =
        (!x_oob(i + 1, j, k) && faceXLiquidWeights.fieldNC()->getValue(i + 1, j, k)) && (!x_oob(i, j, k) && faceXLiquidWeights.fieldNC()->getValue(i, j, k)) &&
        (!y_oob(i, j + 1, k) && faceYLiquidWeights.fieldNC()->getValue(i, j + 1, k)) && (!y_oob(i, j, k) && faceYLiquidWeights.fieldNC()->getValue(i, j, k)) &&
        (!z_oob(i, j, k + 1) && faceZLiquidWeights.fieldNC()->getValue(i, j, k + 1)) && (!z_oob(i, j, k) && faceZLiquidWeights.fieldNC()->getValue(i, j, k));

    if (insystem)
        return MaterialLabels::GENERICFLUID;
    else
        return MaterialLabels::UNSOLVED;
}

void
HDK_PolyStokes::Solver::constructReducedRegions()
{
    // in: centerLabels = MaterialLabels::{SOLID, GENERICFLUID, UNSOLVED}
    // out: centerLabels = MaterialLabels::{SOLID, ACTIVEFLUID, REDUCED, UNSOLVED}

    constructAirBoundaryLayer();
    constructSolidBoundaryLayer();
    if (myDoTile)
        constructTiles();
    overwriteIndices(centerLabels, MaterialLabels::GENERICFLUID, MaterialLabels::REDUCED);
}

void
HDK_PolyStokes::Solver::constructOnlyActiveRegions()
{
    // in: centerLabels = MaterialLabels::{SOLID, GENERICFLUID, UNSOLVED}
    // out: centerLabels = MaterialLabels::{SOLID, ACTIVEFLUID, UNSOLVED}

    overwriteIndices(centerLabels, MaterialLabels::GENERICFLUID, MaterialLabels::ACTIVEFLUID);
}

void
HDK_PolyStokes::Solver::classifyFaces()
{
    classifyFaceAxis(faceXLabels, SamplingType::FACEX);
    classifyFaceAxis(faceYLabels, SamplingType::FACEY);
    classifyFaceAxis(faceZLabels, SamplingType::FACEZ);
}

void
HDK_PolyStokes::Solver::classifyEdges()
{
    classifyEdgeAxis(edgeXYLabels, SamplingType::EDGEXY);
    classifyEdgeAxis(edgeXZLabels, SamplingType::EDGEXZ);
    classifyEdgeAxis(edgeYZLabels, SamplingType::EDGEYZ);
}

void
HDK_PolyStokes::Solver::constructCenterReducedIndices()
{
    const SIM_RawField* faceLiquidWeights[3] = { liquidWeights[4], liquidWeights[5], liquidWeights[6] };

    SIM_VolumetricConnectedComponentBuilder reducedRegionBuilder(centerReducedIndices, centerLabels, faceLiquidWeights);
    myInteriorRegionCount = reducedRegionBuilder.buildConnectedComponents([](const exint label) {
        return label == MaterialLabels::REDUCED;
        });

    overwriteIndices(centerReducedIndices,
        SIM_VolumetricConnectedComponentBuilder::INACTIVE_REGION,
        MaterialLabels::UNASSIGNED);

    // make sure no two stencils overlap
    fixReducedRegionBoundaries();
    // create bounding boxes and eliminate any regions that are too small
    fixSmallReducedRegions();

    // now we are done setting up the indices for the cell centers
    // since we now know the number of reduced regions, create the grain size
    myGrainSize = myInteriorRegionCount / (4 * myThreadCount);
}

void
HDK_PolyStokes::Solver::constructFacesReducedIndices()
{
    constructFaceAxisReducedIndices(faceXReducedIndices, faceXLabels, SamplingType::FACEX);
    constructFaceAxisReducedIndices(faceYReducedIndices, faceYLabels, SamplingType::FACEY);
    constructFaceAxisReducedIndices(faceZReducedIndices, faceZLabels, SamplingType::FACEZ);
}

void
HDK_PolyStokes::Solver::constructEdgesReducedIndices()
{
    constructEdgeAxisReducedIndices(edgeXYReducedIndices, edgeXYLabels, SamplingType::EDGEXY);
    constructEdgeAxisReducedIndices(edgeXZReducedIndices, edgeXZLabels, SamplingType::EDGEXZ);
    constructEdgeAxisReducedIndices(edgeYZReducedIndices, edgeYZLabels, SamplingType::EDGEYZ);
}

void
HDK_PolyStokes::Solver::constructCenterActiveIndices()
{
    overwriteIndices(centerLabels, MaterialLabels::GENERICFLUID, MaterialLabels::ACTIVEFLUID);
    nCenter = serialAssignFieldIndices(centerActiveIndices, centerLabels, 0);
}

void
HDK_PolyStokes::Solver::constructFacesActiveIndices()
{
    overwriteIndices(faceXLabels, MaterialLabels::GENERICFLUID, MaterialLabels::ACTIVEFLUID);
    overwriteIndices(faceYLabels, MaterialLabels::GENERICFLUID, MaterialLabels::ACTIVEFLUID);
    overwriteIndices(faceZLabels, MaterialLabels::GENERICFLUID, MaterialLabels::ACTIVEFLUID);
    nFaceX = serialAssignFieldIndices(faceXActiveIndices, faceXLabels, 0);
    nFaceY = serialAssignFieldIndices(faceYActiveIndices, faceYLabels, 0);
    nFaceZ = serialAssignFieldIndices(faceZActiveIndices, faceZLabels, 0);
}

void
HDK_PolyStokes::Solver::constructEdgesActiveIndices()
{
    overwriteIndices(edgeXYLabels, MaterialLabels::GENERICFLUID, MaterialLabels::ACTIVEFLUID);
    overwriteIndices(edgeXZLabels, MaterialLabels::GENERICFLUID, MaterialLabels::ACTIVEFLUID);
    overwriteIndices(edgeYZLabels, MaterialLabels::GENERICFLUID, MaterialLabels::ACTIVEFLUID);
    nEdgeXY = serialAssignFieldIndices(edgeXYActiveIndices, edgeXYLabels, 0);
    nEdgeXZ = serialAssignFieldIndices(edgeXZActiveIndices, edgeXZLabels, 0);
    nEdgeYZ = serialAssignFieldIndices(edgeYZActiveIndices, edgeYZLabels, 0);
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// called in constructReducedRegions()
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
HDK_PolyStokes::Solver::constructAirBoundaryLayer()
{
    auto cellCompare = [&](const UT_Vector3I& cellA, const UT_Vector3I& cellB)
    {
        // Compare tile number first
        int tileNumberA = centerLabels.field()->indexToLinearTile(cellA[0], cellA[1], cellA[2]);
        int tileNumberB = centerLabels.field()->indexToLinearTile(cellB[0], cellB[1], cellB[2]);

        if (tileNumberA < tileNumberB)
            return true;
        else if (tileNumberA == tileNumberB)
        {
            if (cellA[2] < cellB[2])
                return true;
            else if (cellA[2] == cellB[2])
            {
                if (cellA[1] < cellB[1])
                    return true;
                else if (cellA[1] == cellB[1] &&
                    cellA[0] < cellB[0])
                    return true;
            }
        }

        return false;
    };

    UT_Array<UT_Array<UT_Vector3I>> parallelActiveCenterList;
    parallelActiveCenterList.setSize(myThreadCount);
    buildInitialAirBoundaryLayer(parallelActiveCenterList);
    UT_Array<UT_Vector3I> serialActiveCenterList;

    UT_Array<bool> isTileOccupiedList;
    isTileOccupiedList.setSize(centerLabels.field()->numTiles());

    // Now flood inwards one layer at a time, repeating the process
    for (int layer = 0; layer < myLiquidBoundaryLayerSize - 1; ++layer)
    {
        // Combine parallel active cell lists
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelActiveCenterList[thread].size();

        serialActiveCenterList.clear();
        serialActiveCenterList.bumpCapacity(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
        {
            serialActiveCenterList.concat(parallelActiveCenterList[thread]);
            parallelActiveCenterList[thread].clear();
        }

        UTparallelSort(serialActiveCenterList.begin(), serialActiveCenterList.end(), cellCompare);

        // Build tile list
        isTileOccupiedList.constant(false);
        findOccupiedIndexTiles(isTileOccupiedList, serialActiveCenterList, centerLabels);
        // Uncompress tiles of the new active cells
        uncompressTiles(centerLabels, isTileOccupiedList);

        setActiveLayerCells(serialActiveCenterList);

        // Build next layer of active cells unless we're at the final layer
        if (layer < myLiquidBoundaryLayerSize - 2)
        {
            // TODO: sort list and prevent duplicates
            buildNextLiquidBoundaryLayer(parallelActiveCenterList, serialActiveCenterList);
        }
    }

}

void
HDK_PolyStokes::Solver::buildInitialAirBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>>& parallelActiveCenterList)
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildInitialActiveCellLayerAlgorithm;
    buildInitialActiveCellLayerAlgorithm.run([&](const UT_JobInfo& info)
        {
            UT_Array<UT_Vector3I>& localActiveCenterList = parallelActiveCenterList[info.job()];

            UT_VoxelArrayIteratorI vit;
            vit.setConstArray(centerLabels.field());
            vit.splitByTile(info);

            UT_VoxelTileIteratorI vitt;

            for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
            {
                if (boss->opInterrupt())
                    break;

                if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::GENERICFLUID)
                {
                    vitt.setTile(vit);

                    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                    {
                        if (vitt.getValue() == MaterialLabels::GENERICFLUID)
                        {
                            UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

                            bool isBoundaryCell = false;
                            if (!isBoundaryCell)
                                for (int axis : {0, 1, 2})
                                    for (int direction : {0, 1})
                                    {
                                        UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

                                        if (adjacentCell[axis] < 0 || adjacentCell[axis] >= centerLabels.getVoxelRes()[axis])
                                            continue;

                                        UT_Vector3I face = cellToFaceMap(cell, axis, direction);

                                        // we already know this cell is in the solver, so we can just check if it's bordering an
                                            // unsolved cell and we'll know it's on the boundary
                                        if (getFieldValue(centerLabels, adjacentCell) == MaterialLabels::UNSOLVED)
                                            isBoundaryCell = true;
                                        // if one of the neighbouring faces is not fully fluid, then we also know that this cell
                                            // is on the boundary
                                        if (getFieldValue(*liquidWeights[(uint)faceAxisToSamplingType(axis)], face) < 1.)
                                            isBoundaryCell = true;
                                    }

                            if (isBoundaryCell)
                                localActiveCenterList.append(cell);
                        }
                    }
                }
            }

            return 0;
        });
}

void
HDK_PolyStokes::Solver::buildNextLiquidBoundaryLayer(
    UT_Array<UT_Array<UT_Vector3I> >& parallelNewActiveCellLayer,
    const UT_Array<UT_Vector3I>& oldActiveCellLayer)
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildNextActiveCellLayerAlgorithm;
    buildNextActiveCellLayerAlgorithm.run([&](const UT_JobInfo& info)
        {
            if (boss->opInterrupt())
                return 0;

            UT_Array<UT_Vector3I>& localNewActiveCellLayer = parallelNewActiveCellLayer[info.job()];

            exint startIndex, endIndex;
            const exint elementSize = oldActiveCellLayer.entries();
            info.divideWork(elementSize, startIndex, endIndex);

            if (startIndex > 0)
            {
                while (startIndex != endIndex && oldActiveCellLayer[startIndex] == oldActiveCellLayer[startIndex - 1])
                    ++startIndex;
            }

            UT_Vector3I oldCell(-1, -1, -1);

            const exint localEndIndex = endIndex;
            for (exint cellIndex = startIndex; cellIndex < localEndIndex; ++cellIndex)
            {
                UT_Vector3I cell = oldActiveCellLayer[cellIndex];

                if (cell == oldCell)
                    continue;

                oldCell = cell;

                //assert(isInSystem(myCentralIndex(cell.x(), cell.y(), cell.z())));

                for (int axis : {0, 1, 2}) {
                    // oob checks
                    int n;
                    switch (axis)
                    {
                    case 0:
                        n = nx; break;
                    case 1:
                        n = ny; break;
                    case 2:
                        n = nz; break;
                    }

                    for (int direction : {0, 1})
                    {
                        UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

                        if (adjacentCell[axis] < 0 || adjacentCell[axis] >= n)
                            continue;

                        UT_Vector3I face = cellToFaceMap(cell, axis, direction);

                        const UT_VoxelArrayF& face_weights = *liquidWeights[(int)faceAxisToSamplingType(axis)]->field();

                        if (face_weights(face.x(), face.y(), face.z()) > 0. &&
                            getFieldValue(centerLabels, adjacentCell) == MaterialLabels::GENERICFLUID)
                            localNewActiveCellLayer.append(adjacentCell);
                    }
                }
            }

            return 0;
        });
}

void
HDK_PolyStokes::Solver::constructSolidBoundaryLayer()
{
    UT_Array<UT_Array<UT_Vector3I>> parallelActiveCenterList;
    parallelActiveCenterList.setSize(myThreadCount);
    buildInitialSolidBoundaryLayer(parallelActiveCenterList);
    UT_Array<UT_Vector3I> serialActiveCenterList;

    SIM_RawIndexField visitedCells;
    visitedCells.match(centerLabels);
    visitedCells.makeConstant(MaterialLabels::UNVISITED);

    UT_Array<bool> isTileOccupiedList;
    isTileOccupiedList.setSize(centerLabels.field()->numTiles());

    // Now flood inwards one layer at a time, repeating the process
    for (int layer = 0; layer < mySolidBoundaryLayerSize; ++layer)
    {
        // Combine parallel active cell lists
        exint listSize = 0;
        for (int thread = 0; thread < myThreadCount; ++thread)
            listSize += parallelActiveCenterList[thread].size();

        serialActiveCenterList.clear();
        serialActiveCenterList.bumpCapacity(listSize);

        for (int thread = 0; thread < myThreadCount; ++thread)
        {
            serialActiveCenterList.concat(parallelActiveCenterList[thread]);
            parallelActiveCenterList[thread].clear();
        }

        // Build tile list
        isTileOccupiedList.constant(false);
        findOccupiedIndexTiles(isTileOccupiedList, serialActiveCenterList, centerLabels);

        // Uncompress tiles of the new active cells
        uncompressTiles(centerLabels, isTileOccupiedList);
        uncompressTiles(visitedCells, isTileOccupiedList);

        // Assign active cells
        UTparallelForLightItems(UT_BlockedRange<exint>(0, serialActiveCenterList.size()), [&](const UT_BlockedRange<exint>& range)
            {
                using SIM::FieldUtils::setFieldValue;
                for (exint i = range.begin(); i != range.end(); ++i)
                {
                    UT_Vector3I cell = serialActiveCenterList[i];
                    setFieldValue(centerLabels, cell, MaterialLabels::ACTIVEFLUID);
                    setFieldValue(visitedCells, cell, MaterialLabels::VISITED);
                }
            });

        // Build next layer of active cells unless we're at the final layer
        if (layer < mySolidBoundaryLayerSize - 1)
        {
            // TODO: sort list and prevent duplicates
            buildNextSolidBoundaryLayer(parallelActiveCenterList,
                serialActiveCenterList,
                visitedCells);
        }
    }
}

void
HDK_PolyStokes::Solver::buildInitialSolidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>>& parallelActiveCenterLis)
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildInitialActiveCellLayerAlgorithm;
    buildInitialActiveCellLayerAlgorithm.run([&](const UT_JobInfo& info)
        {
            UT_Array<UT_Vector3I>& localActiveCellList = parallelActiveCenterLis[info.job()];

            UT_VoxelArrayIteratorI vit;
            vit.setConstArray(centerLabels.field());
            vit.splitByTile(info);

            UT_VoxelTileIteratorI vitt;

            for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
            {
                if (boss->opInterrupt())
                    break;

                if (!vit.isTileConstant() ||
                    vit.getValue() == MaterialLabels::GENERICFLUID ||
                    vit.getValue() == MaterialLabels::ACTIVEFLUID)
                {
                    vitt.setTile(vit);

                    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                    {
                        if (vitt.getValue() == MaterialLabels::GENERICFLUID ||
                            vitt.getValue() == MaterialLabels::ACTIVEFLUID)
                        {
                            UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

                            bool isBoundaryCell = false;
                            for (int axis : {0, 1, 2})
                                for (int direction : {0, 1})
                                {
                                    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

                                    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= centerLabels.getVoxelRes()[axis])
                                    {
                                        isBoundaryCell = true;
                                        continue;
                                    }

                                    UT_Vector3I face = cellToFaceMap(cell, axis, direction);

                                    // here we already know this cell is int he solver, so just check if it's
                                        // next to a solid cell
                                    if (//getFieldValue(faceVolumes[axis], face) < 1 &&
                                        getFieldValue(centerLabels, adjacentCell) == MaterialLabels::SOLID)
                                        isBoundaryCell = true;
                                }

                            if (isBoundaryCell)
                                localActiveCellList.append(cell);
                        }
                    }
                }
            }

            return 0;
        });
}

void
HDK_PolyStokes::Solver::buildNextSolidBoundaryLayer(
    UT_Array<UT_Array<UT_Vector3I> >& parallelNewActiveCellLayer,
    const UT_Array<UT_Vector3I>& oldActiveCellLayer,
    const SIM_RawIndexField& visitedCells)
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildNextActiveCellLayerAlgorithm;
    buildNextActiveCellLayerAlgorithm.run([&](const UT_JobInfo& info)
        {
            UT_Array<UT_Vector3I>& localActiveCellList = parallelNewActiveCellLayer[info.job()];

            exint start, end;
            const exint elementSize = oldActiveCellLayer.entries();
            info.divideWork(elementSize, start, end);

            if (boss->opInterrupt())
                return 0;

            const exint localEnd = end;
            for (exint i = start; i < localEnd; ++i)
            {
                if (!(i & 127))
                {
                    if (boss->opInterrupt())
                        return 0;
                }

                UT_Vector3I cell = oldActiveCellLayer[i];

                assert(getFieldValue(centerLabels, cell) == MaterialLabels::ACTIVEFLUID);
                assert(getFieldValue(visitedCells, cell) == MaterialLabels::VISITED);

                for (int axis : {0, 1, 2})
                    for (int direction : {0, 1})
                    {
                        UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

                        if (adjacentCell[axis] < 0 || adjacentCell[axis] >= centerLabels.getVoxelRes()[axis])
                            continue;

                        UT_Vector3I face = cellToFaceMap(cell, axis, direction);

                        if (getFieldValue(*liquidWeights[(int)faceAxisToSamplingType(axis)], face) > 0)
                        {
                            if (getFieldValue(visitedCells, adjacentCell) == MaterialLabels::UNVISITED &&
                                (getFieldValue(centerLabels, adjacentCell) == MaterialLabels::ACTIVEFLUID ||
                                    getFieldValue(centerLabels, adjacentCell) == MaterialLabels::GENERICFLUID))
                                localActiveCellList.append(adjacentCell);
                        }
                    }
            }

            return 0;
        });
}

void
HDK_PolyStokes::Solver::constructTiles()
{
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();

    UTparallelForEachNumber(centerLabels.field()->numTiles(), [&](const UT_BlockedRange<int>& range)
        {
            UT_VoxelArrayIteratorI vit;
            vit.setConstArray(centerLabels.field());

            if (boss->opInterrupt())
                return;

            for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
            {
                vit.myTileStart = tileNumber;
                vit.myTileEnd = tileNumber + 1;
                vit.rewind();

                if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::GENERICFLUID)
                {
                    for (; !vit.atEnd(); vit.advance())
                    {
                        if (vit.getValue() == MaterialLabels::GENERICFLUID)
                        {
                            UT_Vector3I cell(vit.x(), vit.y(), vit.z());

                            for (int i = 0; i < myTilePadding; i++)
                            {
                                if (cell[0] % myTileSize == i ||
                                    cell[1] % myTileSize == i ||
                                    cell[2] % myTileSize == i)
                                    setFieldValue(centerLabels, cell, MaterialLabels::ACTIVEFLUID);
                            }
                        }
                    }
                }
            }
        });
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// called in classifyFaces()
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
HDK_PolyStokes::Solver::classifyFaceAxisPartial(
    SIM_RawIndexField& faceLabels,
    SamplingType type,
    const UT_JobInfo& info
)
{
    using SIM::FieldUtils::getFieldValue;

    UT_VoxelArrayIteratorI vit(faceLabels.fieldNC());
    vit.splitByTile(info);

    for (vit.rewind(); !vit.atEnd(); vit.advance())
    {
        int i = vit.x(), j = vit.y(), k = vit.z();
        UT_Vector3I face(vit.x(), vit.y(), vit.z());

        MaterialLabels label = MaterialLabels::UNASSIGNED;
        label = findFaceLabelFromCenter(i, j, k, type);
        //label = findFaceLabelFromCenterAlt(i, j, k, type);
        /*
        if (getFieldValue(*faceFluidWeights(samplingTypeToFaceAxis(type)), face) > 0.
            && getFieldValue(*faceLiquidWeights(samplingTypeToFaceAxis(type)), face) > 0.)
            label = findFaceLabelFromCenterAlt(i, j, k, type);
            //label = findFaceLabelFromCenter(i, j, k, type);
        else
            label = MaterialLabels::UNSOLVED;
        */
        faceLabels.fieldNC()->setValue(UT_Vector3I(i, j, k), label);
    }
}

HDK_PolyStokes::Solver::MaterialLabels
HDK_PolyStokes::Solver::findFaceLabelFromCenter(uint i, uint j, uint k, SamplingType type)
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::faceToCellMap;
    using SIM::FieldUtils::faceToEdgeMap;

    MaterialLabels retval = MaterialLabels::UNSOLVED;
    UT_Vector3I face(i, j, k);
    uint axis = samplingTypeToFaceAxis(type);

    bool isActiveVelocity = false;

    for (int direction : {0, 1})
    {
        UT_Vector3i cell = faceToCellMap(face, axis, direction);
        if (c_oob(cell.x(), cell.y(), cell.z()))
            continue;
        if (getFieldValue(centerLiquidWeights, cell) > 0.)
            isActiveVelocity = true;
    }

    if (!isActiveVelocity)
        for (int edgeAxis = 0; edgeAxis < 3 && !isActiveVelocity; ++edgeAxis)
        {
            if (edgeAxis == axis) continue;

            for (int direction : {0, 1})
            {
                UT_Vector3i edge = faceToEdgeMap(face, axis, edgeAxis, direction);

                if (getFieldValue(*edgeLiquidWeights(edgeAxis), edge) > 0.)
                {
                    isActiveVelocity = true;
                    break;
                }
            }
        }

    if (isActiveVelocity)
    {
        if (getFieldValue(*faceFluidWeights(axis), face) < 0.5)
            retval = MaterialLabels::SOLID;
        else
            retval = MaterialLabels::GENERICFLUID;
    }

    return retval;
}

HDK_PolyStokes::Solver::MaterialLabels
HDK_PolyStokes::Solver::findFaceLabelFromCenterAlt(uint i, uint j, uint k, SamplingType type)
{
    bool insystem = false;
    // todo convert fieldNC (nonconst) to field
    switch (type)
    {
    case SamplingType::FACEX:
        insystem = (!x_oob(i, j, k) && faceXFluidWeights.fieldNC()->getValue(i, j, k) && faceXLiquidWeights.fieldNC()->getValue(i, j, k));
        break;
    case SamplingType::FACEY:
        insystem = (!y_oob(i, j, k) && faceYFluidWeights.fieldNC()->getValue(i, j, k) && faceYLiquidWeights.fieldNC()->getValue(i, j, k));
        break;
    case SamplingType::FACEZ:
        insystem = (!z_oob(i, j, k) && faceZFluidWeights.fieldNC()->getValue(i, j, k) && faceZLiquidWeights.fieldNC()->getValue(i, j, k));
        break;
    }

    if (!insystem)
        return MaterialLabels::UNSOLVED;

    // check solid boundaries
    switch (type)
    {
    case SamplingType::FACEX:
        if (x_oob(i + 1, j, k) || x_oob(i - 1, j, k))
            return MaterialLabels::SOLID;
        if (faceXFluidWeights.fieldNC()->getValue(i, j, k) < 0.5 || c_oob(i, j, k) || !centerFluidWeights.fieldNC()->getValue(i, j, k) || c_oob(i - 1, j, k) || !centerFluidWeights.fieldNC()->getValue(i - 1, j, k) ||
            !edgeXZFluidWeights.fieldNC()->getValue(i, j, k) || !edgeXZFluidWeights.fieldNC()->getValue(i, j, k + 1) ||
            !edgeXYFluidWeights.fieldNC()->getValue(i, j, k) || !edgeXYFluidWeights.fieldNC()->getValue(i, j + 1, k))
            return MaterialLabels::SOLID;
        break;
    case SamplingType::FACEY:
        if (y_oob(i, j + 1, k) || y_oob(i, j - 1, k))
            return MaterialLabels::SOLID;
        if (faceYFluidWeights.fieldNC()->getValue(i, j, k) < 0.5 || c_oob(i, j, k) || !centerFluidWeights.fieldNC()->getValue(i, j, k) || c_oob(i, j - 1, k) || !centerFluidWeights.fieldNC()->getValue(i, j - 1, k) ||
            !edgeYZFluidWeights.fieldNC()->getValue(i, j, k) || !edgeYZFluidWeights.fieldNC()->getValue(i, j, k + 1) ||
            !edgeXYFluidWeights.fieldNC()->getValue(i, j, k) || !edgeXYFluidWeights.fieldNC()->getValue(i + 1, j, k))
            return MaterialLabels::SOLID;
        break;
    case SamplingType::FACEZ:
        if (z_oob(i, j, k + 1) || z_oob(i, j, k - 1))
            return MaterialLabels::SOLID;
        if (faceZFluidWeights.fieldNC()->getValue(i, j, k) < 0.5 || c_oob(i, j, k) || !centerFluidWeights.fieldNC()->getValue(i, j, k) || c_oob(i, j, k - 1) || !centerFluidWeights.fieldNC()->getValue(i, j, k - 1) ||
            !edgeYZFluidWeights.fieldNC()->getValue(i, j, k) || !edgeYZFluidWeights.fieldNC()->getValue(i, j + 1, k) ||
            !edgeXZFluidWeights.fieldNC()->getValue(i, j, k) || !edgeXZFluidWeights.fieldNC()->getValue(i + 1, j, k))
            return MaterialLabels::SOLID;
        break;
    }

    switch (type)
    {
    case SamplingType::FACEX:
        insystem =
            (centerFluidWeights.fieldNC()->getValue(i, j, k) &&
                edgeXYFluidWeights.fieldNC()->getValue(i, j, k) && !txy_oob(i, j + 1, k) && edgeXYFluidWeights.fieldNC()->getValue(i, j + 1, k) &&
                edgeXZFluidWeights.fieldNC()->getValue(i, j, k) && !txz_oob(i, j, k + 1) && edgeXZFluidWeights.fieldNC()->getValue(i, j, k + 1));
        break;

    case SamplingType::FACEY:
        insystem =
            (centerFluidWeights.fieldNC()->getValue(i, j, k) &&
                edgeXYFluidWeights.fieldNC()->getValue(i, j, k) && !txy_oob(i + 1, j, k) && edgeXYFluidWeights.fieldNC()->getValue(i + 1, j, k) &&
                edgeYZFluidWeights.fieldNC()->getValue(i, j, k) && !tyz_oob(i, j, k + 1) && edgeYZFluidWeights.fieldNC()->getValue(i, j, k + 1));
        break;

    case SamplingType::FACEZ:
        insystem =
            (centerFluidWeights.fieldNC()->getValue(i, j, k) &&
                edgeXZFluidWeights.fieldNC()->getValue(i, j, k) && !txz_oob(i + 1, j, k) && edgeXZFluidWeights.fieldNC()->getValue(i + 1, j, k) &&
                edgeYZFluidWeights.fieldNC()->getValue(i, j, k) && !tyz_oob(i, j + 1, k) && edgeYZFluidWeights.fieldNC()->getValue(i, j + 1, k));
        break;
    }

    if (insystem)
        return MaterialLabels::GENERICFLUID;
    else
        return MaterialLabels::UNSOLVED;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// called in classifyEdges()
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
HDK_PolyStokes::Solver::classifyEdgeAxisPartial(
    SIM_RawIndexField& edgeLabels,
    SamplingType type,
    const UT_JobInfo& info
)
{
    using SIM::FieldUtils::getFieldValue;

    UT_VoxelArrayIteratorI vit(edgeLabels.fieldNC());
    vit.splitByTile(info);

    for (vit.rewind(); !vit.atEnd(); vit.advance())
    {
        int i = vit.x(), j = vit.y(), k = vit.z();
        UT_Vector3I edge(vit.x(), vit.y(), vit.z());

        MaterialLabels label = MaterialLabels::UNASSIGNED;
        // in progress clamping fluid weight
        //if (getFieldValue(*edgeFluidWeights(samplingTypeToEdgeAxis(type)), edge) > 0.
            //&& getFieldValue(*edgeLiquidWeights(samplingTypeToEdgeAxis(type)), edge) > 0.)

        /*if (getFieldValue(*edgeLiquidWeights(samplingTypeToEdgeAxis(type)), edge) > 0.)
            label = findEdgeLabelFromFace(i, j, k, type);
        else
            label = MaterialLabels::UNSOLVED;
        */
        label = findEdgeLabelFromFaceAlt(i, j, k, type);
        edgeLabels.fieldNC()->setValue(UT_Vector3I(i, j, k), label);
    }
}

HDK_PolyStokes::Solver::MaterialLabels
HDK_PolyStokes::Solver::findEdgeLabelFromFace(uint i, uint j, uint k, SamplingType type)
{
    MaterialLabels retval = MaterialLabels::UNSOLVED;

    switch (type)
    {
    case SamplingType::EDGEXY:
        // check if there should be an active DOF here
        if (!x_oob(i, j, k) && isSolved(faceXLabels(i, j, k)))
            retval = MaterialLabels::GENERICFLUID;
        else if (!x_oob(i, j - 1, k) && isSolved(faceXLabels(i, j - 1, k)))
            retval = MaterialLabels::GENERICFLUID;
        else if (!y_oob(i, j, k) && isSolved(faceYLabels(i, j, k)))
            retval = MaterialLabels::GENERICFLUID;
        else if (!y_oob(i - 1, j, k) && isSolved(faceYLabels(i - 1, j, k)))
            retval = MaterialLabels::GENERICFLUID;

        // in progress clamping for fluid weight
        // check if this should be a solid BC
        if (retval == MaterialLabels::GENERICFLUID && edgeXYFluidWeights.fieldNC()->getValue(i, j, k) == 0.)
            retval = MaterialLabels::SOLID;
        if (txy_oob(i - 1, j, k) || txy_oob(i + 1, j, k) || txy_oob(i, j - 1, k) || txy_oob(i, j + 1, k))
            retval = MaterialLabels::SOLID;

        break;
    case SamplingType::EDGEXZ:
        // check if there should be an active DOF here
        if (!x_oob(i, j, k) && isSolved(faceXLabels(i, j, k)))
            retval = MaterialLabels::GENERICFLUID;
        else if (!x_oob(i, j, k - 1) && isSolved(faceXLabels(i, j, k - 1)))
            retval = MaterialLabels::GENERICFLUID;
        else if (!z_oob(i, j, k) && isSolved(faceZLabels(i, j, k)))
            retval = MaterialLabels::GENERICFLUID;
        else if (!z_oob(i - 1, j, k) && isSolved(faceZLabels(i - 1, j, k)))
            retval = MaterialLabels::GENERICFLUID;

        // check if this should be a solid BC
        if (retval == MaterialLabels::GENERICFLUID && edgeXZFluidWeights.fieldNC()->getValue(i, j, k) == 0.)
            retval = MaterialLabels::SOLID;
        if (txz_oob(i - 1, j, k) || txz_oob(i + 1, j, k) || txz_oob(i, j, k - 1) || txz_oob(i, j, k + 1))
            retval = MaterialLabels::SOLID;

        break;
    case SamplingType::EDGEYZ:
        // check if there should be an active DOF here
        if (!y_oob(i, j, k) && isSolved(faceYLabels(i, j, k)))
            retval = MaterialLabels::GENERICFLUID;
        else if (!y_oob(i, j, k - 1) && isSolved(faceYLabels(i, j, k - 1)))
            retval = MaterialLabels::GENERICFLUID;
        else if (!z_oob(i, j, k) && isSolved(faceZLabels(i, j, k)))
            retval = MaterialLabels::GENERICFLUID;
        else if (!z_oob(i, j - 1, k) && isSolved(faceZLabels(i, j - 1, k)))
            retval = MaterialLabels::GENERICFLUID;

        // check if this should be a solid BC
        if (retval == MaterialLabels::GENERICFLUID && edgeYZFluidWeights.fieldNC()->getValue(i, j, k) == 0.)
            retval = MaterialLabels::SOLID;
        if (tyz_oob(i, j - 1, k) || tyz_oob(i, j + 1, k) || tyz_oob(i, j, k - 1) || tyz_oob(i, j, k + 1))
            retval = MaterialLabels::SOLID;

        break;
    default:
        assert(type == SamplingType::EDGEXY || type == SamplingType::EDGEXZ || type == SamplingType::EDGEYZ);
        break;
    }

    return retval;
}

HDK_PolyStokes::Solver::MaterialLabels
HDK_PolyStokes::Solver::findEdgeLabelFromFaceAlt(uint i, uint j, uint k, SamplingType type)
{
    bool insystem = false;

    switch (type)
    {
    case SamplingType::EDGEXY:
        insystem = (!txy_oob(i, j, k) && edgeXYLiquidWeights.fieldNC()->getValue(i, j, k) && edgeXYFluidWeights.fieldNC()->getValue(i, j, k));
        break;

    case SamplingType::EDGEXZ:
        insystem = (!txz_oob(i, j, k) && edgeXZLiquidWeights.fieldNC()->getValue(i, j, k) && edgeXZFluidWeights.fieldNC()->getValue(i, j, k));
        break;

    case SamplingType::EDGEYZ:
        insystem = (!tyz_oob(i, j, k) && edgeYZLiquidWeights.fieldNC()->getValue(i, j, k) && edgeYZFluidWeights.fieldNC()->getValue(i, j, k));
        break;
    }

    if (!insystem)
        return MaterialLabels::UNSOLVED;

    switch (type)
    {
    case SamplingType::EDGEXY:
        insystem =
            faceXLiquidWeights.fieldNC()->getValue(i, j, k) && !x_oob(i, j - 1, k) && faceXLiquidWeights.fieldNC()->getValue(i, j - 1, k) &&
            faceYLiquidWeights.fieldNC()->getValue(i, j, k) && !y_oob(i - 1, j, k) && faceYLiquidWeights.fieldNC()->getValue(i - 1, j, k);
        break;
    case SamplingType::EDGEXZ:
        insystem =
            faceXLiquidWeights.fieldNC()->getValue(i, j, k) && !x_oob(i, j, k - 1) && faceXLiquidWeights.fieldNC()->getValue(i, j, k - 1) &&
            faceZLiquidWeights.fieldNC()->getValue(i, j, k) && !z_oob(i - 1, j, k) && faceZLiquidWeights.fieldNC()->getValue(i - 1, j, k);
        break;
    case SamplingType::EDGEYZ:
        insystem =
            faceYLiquidWeights.fieldNC()->getValue(i, j, k) && !y_oob(i, j, k - 1) && faceYLiquidWeights.fieldNC()->getValue(i, j, k - 1) &&
            faceZLiquidWeights.fieldNC()->getValue(i, j, k) && !z_oob(i, j - 1, k) && faceZLiquidWeights.fieldNC()->getValue(i, j - 1, k);
        break;
    }

    if (insystem)
        return MaterialLabels::GENERICFLUID;
    else
        return MaterialLabels::UNSOLVED;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// called in constructCenterReducedIndices()
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
HDK_PolyStokes::Solver::fixReducedRegionBoundaries()
{
    // converts any cells which are too close on the boundary between two neighbouring reduced regions to
    // uniform cells, to ensure no two reduced regions' stencils touch

    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(centerLabels.field());
    UT_VoxelTileIteratorI vitt;
    UT_Interrupt* boss = UTgetInterrupt();

    bool done = false;

    int nRemoved = 0;
    int nLoops = 0;

    while (!done)
    {
        done = true;
        nLoops++;

        // todo make sure this iterator works (ie restarts every while loop)
        for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
        {
            if (boss->opInterrupt())
                break;

            if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::ACTIVEFLUID)
            {
                vitt.setTile(vit);

                for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                {
                    if (vitt.getValue() == MaterialLabels::ACTIVEFLUID)
                    {
                        UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());
                        bool applyFix = false;
                        bool isBoundaryCell = false;
                        int adjacentInteriorRegion;

                        // todo exit out of loop once we know we need to apply fix
                        for (int axis : {0, 1, 2})
                            for (int direction : {0, 1})
                            {
                                UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);
                                if (isReduced(getFieldValue(centerLabels, adjacentCell)))
                                    if (isBoundaryCell == false) {
                                        isBoundaryCell = true;
                                        adjacentInteriorRegion = getFieldValue(centerReducedIndices, adjacentCell);
                                    }
                                    else if (getFieldValue(centerReducedIndices, adjacentCell) != adjacentInteriorRegion)
                                        applyFix = true;

                                // check one more star layer in
                                /*
                                for (int axis2 : {0, 1, 2})
                                    for (int direction2 : {0, 1})
                                    {
                                        // or just check corners with this line
                                        if (direction2 == direction) continue;

                                        UT_Vector3I adjacentCell2 = cellToCellMap(adjacentCell, axis2, direction2);
                                        if (adjacentCell2 == cell) continue;
                                        if (isReduced(getFieldValue(centerLabels, adjacentCell2)))
                                            if (isBoundaryCell == false) {
                                                isBoundaryCell = true;
                                                adjacentInteriorRegion = getFieldValue(centerReducedIndices, adjacentCell2);
                                            }
                                            else if (getFieldValue(centerReducedIndices, adjacentCell2) != adjacentInteriorRegion)
                                                applyFix = true;
                                    }
                                */
                            }

                        if (applyFix)
                        {
                            done = false;
                            for (int axis : {0, 1, 2})
                                for (int direction : {0, 1})
                                {
                                    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);
                                    if (isReduced(getFieldValue(centerLabels, adjacentCell)))
                                    {
                                        nRemoved++;
                                        setFieldValue(centerLabels, adjacentCell, MaterialLabels::ACTIVEFLUID);
                                        setFieldValue(centerReducedIndices, adjacentCell, MaterialLabels::UNASSIGNED);
                                    }
                                }
                        }

                    }
                }
            }
        }
    }
}

void
HDK_PolyStokes::Solver::fixSmallReducedRegions()
{
    UT_Array<UT_Array<UT_Vector3I> > parallelInteriorRegionBBMin;
    UT_Array<UT_Array<UT_Vector3I> > parallelInteriorRegionBBMax;

    parallelInteriorRegionBBMin.setSize(myThreadCount);
    parallelInteriorRegionBBMax.setSize(myThreadCount);

    for (int thread = 0; thread < myThreadCount; ++thread)
    {
        parallelInteriorRegionBBMin[thread].setSize(myInteriorRegionCount);
        parallelInteriorRegionBBMin[thread].constant(UT_Vector3I(std::numeric_limits<exint>::max()));

        parallelInteriorRegionBBMax[thread].setSize(myInteriorRegionCount);
        parallelInteriorRegionBBMax[thread].constant(UT_Vector3I(std::numeric_limits<exint>::min()));
    }

    buildInteriorBoundingBoxes(parallelInteriorRegionBBMin, parallelInteriorRegionBBMax);

    UT_Array<UT_Vector3I> interiorRegionBBMin;
    interiorRegionBBMin.setSize(myInteriorRegionCount);
    interiorRegionBBMin.constant(UT_Vector3I(std::numeric_limits<exint>::max()));

    UT_Array<UT_Vector3I> interiorRegionBBMax;
    interiorRegionBBMax.setSize(myInteriorRegionCount);
    interiorRegionBBMax.constant(UT_Vector3I(std::numeric_limits<exint>::min()));

    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
            {
                for (int thread = 0; thread < myThreadCount; ++thread)
                {
                    for (int axis : {0, 1, 2})
                    {
                        interiorRegionBBMin[interiorRegion][axis] = std::min(parallelInteriorRegionBBMin[thread][interiorRegion][axis], interiorRegionBBMin[interiorRegion][axis]);
                        interiorRegionBBMax[interiorRegion][axis] = std::max(parallelInteriorRegionBBMax[thread][interiorRegion][axis], interiorRegionBBMax[interiorRegion][axis]);
                    }
                }
            }
        });

    UT_Array<bool> doRemoveRegion;
    doRemoveRegion.setSize(myInteriorRegionCount);
    doRemoveRegion.constant(false);

    /*
    // shrink to a uniform box here?
    {
        shrinkToBoxReducedRegions(
            interiorRegionBBMin,
            interiorRegionBBMax
        );
    }*/

    tbb::parallel_for(tbb::blocked_range<exint>(0, myInteriorRegionCount, myGrainSize), [&](const tbb::blocked_range<exint>& range)
        {
            for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
            {
                for (int axis : {0, 1, 2})
                {
                    if (interiorRegionBBMax[interiorRegion][axis] == interiorRegionBBMin[interiorRegion][axis])
                        doRemoveRegion[interiorRegion] = true;

                    if (interiorRegionBBMin[interiorRegion][axis] > (interiorRegionBBMax[interiorRegion][axis]-3))
                        doRemoveRegion[interiorRegion] = true;
                }
            }
        });

    UT_Array<exint> remapRegion;
    remapRegion.setSize(myInteriorRegionCount);
    remapRegion.constant(-1);

    int regionMap = 0;
    for (exint interiorRegion = 0; interiorRegion < myInteriorRegionCount; ++interiorRegion)
    {
        if (!doRemoveRegion[interiorRegion])
            remapRegion[interiorRegion] = regionMap++;
    }

    if (myInteriorRegionCount > regionMap)
    {
        myInteriorRegionCount = regionMap;

        remapInteriorRegions(doRemoveRegion, remapRegion);
    }
}

void
HDK_PolyStokes::Solver::remapInteriorRegions(
    const UT_Array<bool>& doRemoveRegion,
    const UT_Array<exint>& remapRegion)
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();

    UTparallelForEachNumber(centerLabels.field()->numTiles(), [&](const UT_BlockedRange<int>& range)
        {
            UT_VoxelArrayIteratorI vit;
            vit.setConstArray(centerLabels.field());

            if (boss->opInterrupt())
                return;

            for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
            {
                vit.myTileStart = tileNumber;
                vit.myTileEnd = tileNumber + 1;
                vit.rewind();

                if (!vit.isTileConstant() || isReduced(vit.getValue()))
                {
                    for (; !vit.atEnd(); vit.advance())
                    {
                        if (isReduced(vit.getValue()))
                        {
                            UT_Vector3I cell(vit.x(), vit.y(), vit.z());

                            exint interiorRegion = getFieldValue(centerReducedIndices, cell);
                            assert(interiorRegion >= 0);

                            if (doRemoveRegion[interiorRegion])
                            {
                                assert(remapRegion[interiorRegion] == -1);
                                setFieldValue(centerLabels, cell, MaterialLabels::ACTIVEFLUID);
                                setFieldValue(centerReducedIndices, cell, MaterialLabels::UNASSIGNED);
                            }
                            else if (remapRegion[interiorRegion] < interiorRegion)
                                setFieldValue(centerReducedIndices, cell, remapRegion[interiorRegion]);
                            else assert(remapRegion[interiorRegion] == interiorRegion);
                        }
                    }
                }
            }
        });
}

void
HDK_PolyStokes::Solver::shrinkToBoxReducedRegions(
    UT_Array<UT_Vector3I>& interiorRegionBBMin,
    UT_Array<UT_Vector3I>& interiorRegionBBMax
)
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();


    for (int interiorRegion = 0; interiorRegion < myInteriorRegionCount; ++interiorRegion)
    {
        bool isAABB = false;
        
        while (!isAABB)
        {
            Index AABBLimits[6];
            for (int axis : {0, 1, 2}) {
                AABBLimits[axis] = interiorRegionBBMin[interiorRegion][axis];
                AABBLimits[axis + 3] = interiorRegionBBMax[interiorRegion][axis];
            }

            // count the number of uniform cells on each side of the box
            Index faceUniformArea[6] = { 0, 0, 0, 0, 0, 0 };
            Index faceArea[6] = { 0, 0, 0, 0, 0, 0 };
            for (int direction : {0, 1}) {
                for (int axis : {0, 1, 2}) {

                    // first construct the face limits
                    Index faceLimits[6];
                    for (int i = 0; i < 6; ++i) faceLimits[i] = AABBLimits[i];              // copy the AABBLimits
                    faceLimits[axis + 3 * (!direction)] = AABBLimits[axis + 3 * direction]; // flatten the face that we are currently on

                    // now loop over this face and count the number of uniform cells
                    Index localUniformArea = 0;
                    for (int i = faceLimits[AABBSide::X_MIN]; i <= faceLimits[AABBSide::X_MAX]; ++i) {
                        for (int j = faceLimits[AABBSide::Y_MIN]; j <= faceLimits[AABBSide::Y_MAX]; ++j) {
                            for (int k = faceLimits[AABBSide::Z_MIN]; k <= faceLimits[AABBSide::Z_MAX]; ++k) {
                                UT_Vector3I cell(i, j, k);
                                if (isActive(getFieldValue(centerLabels, cell))) {
                                    ++localUniformArea;
                                }
                            }
                        }
                    }
                    // get the total area of the face
                    Index localArea =
                        (faceLimits[AABBSide::X_MAX] - faceLimits[AABBSide::X_MIN] + 1)
                        * (faceLimits[AABBSide::Y_MAX] - faceLimits[AABBSide::Y_MIN] + 1)
                        * (faceLimits[AABBSide::Z_MAX] - faceLimits[AABBSide::Z_MIN] + 1);

                    faceUniformArea[axis + direction * 3] = localUniformArea;
                    faceArea[axis + direction * 3] = localArea;
                }
            }

            // now choose which face we want to cull, for now just take the one with the most uniform cells
            AABBSide faceToCull = (AABBSide)0;
            Index currentMaxUniformArea = 0;
            for (int i = 0; i < 6; ++i)
                if (faceUniformArea[i] > currentMaxUniformArea) {
                    faceToCull = (AABBSide)i;
                    currentMaxUniformArea = faceUniformArea[i];
                }

            if (currentMaxUniformArea == 0) {
                isAABB = true;
            }

            // now cull this face by manually setting the grid values
            if (!isAABB) {
                int direction = (faceToCull >= 3) ? 1 : 0;
                int axis = faceToCull % 3;

                // first construct the face limits
                Index faceLimits[6];
                for (int i = 0; i < 6; ++i) faceLimits[i] = AABBLimits[i];              // copy the AABBLimits
                faceLimits[axis + 3 * (!direction)] = AABBLimits[axis + 3 * direction]; // flatten the face that we are currently on

                // now loop over this face and count the number of uniform cells
                Index localUniformArea = 0;
                for (int i = faceLimits[AABBSide::X_MIN]; i <= faceLimits[AABBSide::X_MAX]; ++i)
                    for (int j = faceLimits[AABBSide::Y_MIN]; j <= faceLimits[AABBSide::Y_MAX]; ++j)
                        for (int k = faceLimits[AABBSide::Z_MIN]; k <= faceLimits[AABBSide::Z_MAX]; ++k)
                            if (isReduced(getFieldValue(centerLabels, UT_Vector3I(i, j, k)))) {
                                UT_Vector3I cell(i, j, k);
                                setFieldValue(centerLabels, cell, MaterialLabels::ACTIVEFLUID);
                                setFieldValue(centerReducedIndices, cell, MaterialLabels::UNASSIGNED);
                            }

                // we also need to make sure the bounding box gets updated
                if (direction == 0)
                    ++interiorRegionBBMin[interiorRegion][axis];
                else
                    --interiorRegionBBMax[interiorRegion][axis];
            }
        }
    }

}

void
HDK_PolyStokes::Solver::buildInteriorBoundingBoxes(
    UT_Array<UT_Array<UT_Vector3I> >& parallelInteriorRegionBBMin,
    UT_Array<UT_Array<UT_Vector3I> >& parallelInteriorRegionBBMax)
{
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildInteriorBoundingBoxesAlgorithm;
    buildInteriorBoundingBoxesAlgorithm.run([&](const UT_JobInfo& info)
        {
            UT_Array<UT_Vector3I>& localInteriorRegionBBMin = parallelInteriorRegionBBMin[info.job()];
            UT_Array<UT_Vector3I>& localInteriorRegionBBMax = parallelInteriorRegionBBMax[info.job()];

            UT_VoxelArrayIteratorI vit;
            vit.setConstArray(centerLabels.field());
            vit.splitByTile(info);

            UT_VoxelTileIteratorI vitt;

            for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
            {
                if (boss->opInterrupt())
                    break;

                if (!vit.isTileConstant() || isReduced(vit.getValue()))
                {
                    vitt.setTile(vit);

                    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                    {
                        if (isReduced(vitt.getValue()))
                        {
                            UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());
                            exint interiorRegion = getFieldValue(centerReducedIndices, cell);

                            for (int axis : {0, 1, 2})
                            {
                                localInteriorRegionBBMin[interiorRegion][axis] = std::min(localInteriorRegionBBMin[interiorRegion][axis], cell[axis]);
                                localInteriorRegionBBMax[interiorRegion][axis] = std::max(localInteriorRegionBBMax[interiorRegion][axis], cell[axis]);
                            }
                        }
                    }
                }
            }

            return 0;
        });
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// called in constructFacesReducedIndices()
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
HDK_PolyStokes::Solver::constructFaceAxisReducedIndicesPartial(
    SIM_RawIndexField& faceReducedIndexField,
    SIM_RawIndexField& faceLabels,
    SamplingType type,
    const UT_JobInfo& info
)
{
    // set the tile numbers to the interior field
    UT_VoxelArrayIteratorI vit(faceReducedIndexField.fieldNC());
    //vit.setCompressOnExit(true);
    vit.splitByTile(info);

    for (vit.rewind(); !vit.atEnd(); vit.advance())
    {
        int i = vit.x(), j = vit.y(), k = vit.z();

        auto idx = findFaceReducedIndexFromCenter(i, j, k, type);
        if (idx != MaterialLabels::UNASSIGNED) {
            faceLabels.fieldNC()->setValue(UT_Vector3I(i, j, k), MaterialLabels::REDUCED);
            faceReducedIndexField.fieldNC()->setValue(UT_Vector3I(i, j, k), idx);
        }
    }
}

exint
HDK_PolyStokes::Solver::findFaceReducedIndexFromCenter(uint i, uint j, uint k, SamplingType type)
{
    exint retval = MaterialLabels::UNASSIGNED;

    switch (type)
    {
    case SamplingType::FACEX:
        if (!c_oob(i, j, k) && centerLabels(i, j, k) == MaterialLabels::REDUCED)
            retval = centerReducedIndices(i, j, k);
        else if (!c_oob(i - 1, j, k) && centerLabels(i - 1, j, k) == MaterialLabels::REDUCED)
            retval = centerReducedIndices(i - 1, j, k);
        break;
    case SamplingType::FACEY:
        if (!c_oob(i, j, k) && centerLabels(i, j, k) == MaterialLabels::REDUCED)
            retval = centerReducedIndices(i, j, k);
        else if (!c_oob(i, j - 1, k) && centerLabels(i, j - 1, k) == MaterialLabels::REDUCED)
            retval = centerReducedIndices(i, j - 1, k);
        break;
    case SamplingType::FACEZ:
        if (!c_oob(i, j, k) && centerLabels(i, j, k) == MaterialLabels::REDUCED)
            retval = centerReducedIndices(i, j, k);
        else if (!c_oob(i, j, k - 1) && centerLabels(i, j, k - 1) == MaterialLabels::REDUCED)
            retval = centerReducedIndices(i, j, k - 1);
        break;
    default:
        assert(type == SamplingType::FACEX || type == SamplingType::FACEY || type == SamplingType::FACEZ);
        break;
    }
    return retval;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// called in constructEdgesReducedIndices()
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
HDK_PolyStokes::Solver::constructEdgeAxisReducedIndicesPartial(
    SIM_RawIndexField& edgeReducedIndexField,
    SIM_RawIndexField& edgeLabels,
    SamplingType type,
    const UT_JobInfo& info
)
{
    // set the tile numbers to the interior field
    UT_VoxelArrayIteratorI vit(edgeReducedIndexField.fieldNC());
    //vit.setCompressOnExit(true);
    vit.splitByTile(info);

    for (vit.rewind(); !vit.atEnd(); vit.advance())
    {
        int i = vit.x(), j = vit.y(), k = vit.z();

        /*
        auto idx = findEdgeReducedIndexFromFaces(i, j, k, type);
        if (idx != MaterialLabels::UNASSIGNED) {
            edgeLabels.fieldNC()->setValue(UT_Vector3I(i, j, k), MaterialLabels::REDUCED);
            edgeReducedIndexField.fieldNC()->setValue(UT_Vector3I(i, j, k), idx);
        }*/

        MaterialLabels label = MaterialLabels::UNASSIGNED;
        exint idx = MaterialLabels::UNASSIGNED;

        switch (type)
        {
        case SamplingType::EDGEXY:
            if ((!x_oob(i, j, k) && faceXLabels(i, j, k) == MaterialLabels::REDUCED)
                && (!x_oob(i, j - 1, k) && faceXLabels(i, j - 1, k) == MaterialLabels::REDUCED)
                && (!y_oob(i, j, k) && faceYLabels(i, j, k) == MaterialLabels::REDUCED)
                && (!y_oob(i - 1, j, k) && faceYLabels(i - 1, j, k) == MaterialLabels::REDUCED))
            {
                idx = faceXReducedIndices(i, j, k);
                label = MaterialLabels::REDUCED;
            }
            else if (!x_oob(i, j, k) && faceXLabels(i, j, k) == MaterialLabels::REDUCED)
            {
                idx = faceXReducedIndices(i, j, k);
                label = MaterialLabels::BOUNDARY;
            }
            else if (!x_oob(i, j - 1, k) && faceXLabels(i, j - 1, k) == MaterialLabels::REDUCED)
            {
                idx = faceXReducedIndices(i, j - 1, k);
                label = MaterialLabels::BOUNDARY;
            }
            else if (!y_oob(i, j, k) && faceYLabels(i, j, k) == MaterialLabels::REDUCED)
            {
                idx = faceYReducedIndices(i, j, k);
                label = MaterialLabels::BOUNDARY;
            }
            else if (!y_oob(i - 1, j, k) && faceYLabels(i - 1, j, k) == MaterialLabels::REDUCED)
            {
                idx = faceYReducedIndices(i - 1, j, k);
                label = MaterialLabels::BOUNDARY;
            }
            break;
        case SamplingType::EDGEXZ:
            if ((!x_oob(i, j, k) && faceXLabels(i, j, k) == MaterialLabels::REDUCED)
                && (!x_oob(i, j, k - 1) && faceXLabels(i, j, k - 1) == MaterialLabels::REDUCED)
                && (!z_oob(i, j, k) && faceZLabels(i, j, k) == MaterialLabels::REDUCED)
                && (!z_oob(i - 1, j, k) && faceZLabels(i - 1, j, k) == MaterialLabels::REDUCED))
            {
                idx = faceXReducedIndices(i, j, k);
                label = MaterialLabels::REDUCED;
            }
            else if (!x_oob(i, j, k) && faceXLabels(i, j, k) == MaterialLabels::REDUCED)
            {
                idx = faceXReducedIndices(i, j, k);
                label = MaterialLabels::BOUNDARY;
            }
            else if (!x_oob(i, j, k - 1) && faceXLabels(i, j, k - 1) == MaterialLabels::REDUCED)
            {
                idx = faceXReducedIndices(i, j, k - 1);
                label = MaterialLabels::BOUNDARY;
            }
            else if (!z_oob(i, j, k) && faceZLabels(i, j, k) == MaterialLabels::REDUCED)
            {
                idx = faceZReducedIndices(i, j, k);
                label = MaterialLabels::BOUNDARY;
            }
            else if (!z_oob(i - 1, j, k) && faceZLabels(i - 1, j, k) == MaterialLabels::REDUCED)
            {
                idx = faceZReducedIndices(i - 1, j, k);
                label = MaterialLabels::BOUNDARY;
            }
            break;
        case SamplingType::EDGEYZ:
            if ((!y_oob(i, j, k) && faceYLabels(i, j, k) == MaterialLabels::REDUCED)
                && (!y_oob(i, j, k - 1) && faceYLabels(i, j, k - 1) == MaterialLabels::REDUCED)
                && (!z_oob(i, j, k) && faceZLabels(i, j, k) == MaterialLabels::REDUCED)
                && (!z_oob(i, j - 1, k) && faceZLabels(i, j - 1, k) == MaterialLabels::REDUCED))
            {
                idx = faceYReducedIndices(i, j - 1, k);
                label = MaterialLabels::REDUCED;
            }
            else if (!y_oob(i, j, k) && faceYLabels(i, j, k) == MaterialLabels::REDUCED)
            {
                idx = faceYReducedIndices(i, j, k);
                label = MaterialLabels::BOUNDARY;
            }
            else if (!y_oob(i, j, k - 1) && faceYLabels(i, j, k - 1) == MaterialLabels::REDUCED)
            {
                idx = faceYReducedIndices(i, j, k - 1);
                label = MaterialLabels::BOUNDARY;
            }
            else if (!z_oob(i, j, k) && faceZLabels(i, j, k) == MaterialLabels::REDUCED)
            {
                idx = faceZReducedIndices(i, j, k);
                label = MaterialLabels::BOUNDARY;
            }
            else if (!z_oob(i, j - 1, k) && faceZLabels(i, j - 1, k) == MaterialLabels::REDUCED)
            {
                idx = faceZReducedIndices(i, j - 1, k);
                label = MaterialLabels::BOUNDARY;
            }
            break;
        }
        if (idx != MaterialLabels::UNASSIGNED) {
            edgeLabels.fieldNC()->setValue(UT_Vector3I(i, j, k), label);
            edgeReducedIndexField.fieldNC()->setValue(UT_Vector3I(i, j, k), idx);
        }
    }
}

exint
HDK_PolyStokes::Solver::findEdgeReducedIndexFromFaces(uint i, uint j, uint k, SamplingType type)
{
    exint retval = MaterialLabels::UNASSIGNED;

    switch (type)
    {
    case SamplingType::EDGEXY:
        if ((!x_oob(i, j, k) && faceXLabels(i, j, k) == MaterialLabels::REDUCED)
            && (!x_oob(i, j - 1, k) && faceXLabels(i, j - 1, k) == MaterialLabels::REDUCED)
            && (!y_oob(i, j, k) && faceYLabels(i, j, k) == MaterialLabels::REDUCED)
            && (!y_oob(i - 1, j, k) && faceYLabels(i - 1, j, k) == MaterialLabels::REDUCED))
            retval = faceXReducedIndices(i, j, k);
        break;
    case SamplingType::EDGEXZ:
        if ((!x_oob(i, j, k) && faceXLabels(i, j, k) == MaterialLabels::REDUCED)
            && (!x_oob(i, j, k - 1) && faceXLabels(i, j, k - 1) == MaterialLabels::REDUCED)
            && (!z_oob(i, j, k) && faceZLabels(i, j, k) == MaterialLabels::REDUCED)
            && (!z_oob(i - 1, j, k) && faceZLabels(i - 1, j, k) == MaterialLabels::REDUCED))
            retval = faceXReducedIndices(i, j, k);
        break;
    case SamplingType::EDGEYZ:
        if ((!y_oob(i, j, k) && faceYLabels(i, j, k) == MaterialLabels::REDUCED)
            && (!y_oob(i, j, k - 1) && faceYLabels(i, j, k - 1) == MaterialLabels::REDUCED)
            && (!z_oob(i, j, k) && faceZLabels(i, j, k) == MaterialLabels::REDUCED)
            && (!z_oob(i, j - 1, k) && faceZLabels(i, j - 1, k) == MaterialLabels::REDUCED))
            retval = faceYReducedIndices(i, j - 1, k);
        break;
    default:
        assert(type == SamplingType::EDGEXY || type == SamplingType::EDGEXZ || type == SamplingType::EDGEYZ);
        break;
    }
    /*
        switch (type)
        {
        case SamplingType::EDGEXY:
            if (!x_oob(i, j, k) && faceXLabels(i, j, k) == MaterialLabels::REDUCED)
                retval = faceXReducedIndices(i, j, k);
            else if (!x_oob(i, j - 1, k) && faceXLabels(i, j - 1, k) == MaterialLabels::REDUCED)
                retval = faceXReducedIndices(i, j - 1, k);
            else if (!y_oob(i, j, k) && faceYLabels(i, j, k) == MaterialLabels::REDUCED)
                retval = faceYReducedIndices(i, j, k);
            else if (!y_oob(i - 1, j, k) && faceYLabels(i - 1, j, k) == MaterialLabels::REDUCED)
                retval = faceYReducedIndices(i - 1, j, k);
            break;
        case SamplingType::EDGEXZ:
            if (!x_oob(i, j, k) && faceXLabels(i, j, k) == MaterialLabels::REDUCED)
                retval = faceXReducedIndices(i, j, k);
            else if (!x_oob(i, j, k - 1) && faceXLabels(i, j, k - 1) == MaterialLabels::REDUCED)
                retval = faceXReducedIndices(i, j, k - 1);
            else if (!z_oob(i, j, k) && faceZLabels(i, j, k) == MaterialLabels::REDUCED)
                retval = faceZReducedIndices(i, j, k);
            else if (!z_oob(i - 1, j, k) && faceZLabels(i - 1, j, k) == MaterialLabels::REDUCED)
                retval = faceZReducedIndices(i - 1, j, k);
            break;
        case SamplingType::EDGEYZ:
            if (!y_oob(i, j, k) && faceYLabels(i, j, k) == MaterialLabels::REDUCED)
                retval = faceYReducedIndices(i, j, k);
            else if (!y_oob(i, j, k - 1) && faceYLabels(i, j, k - 1) == MaterialLabels::REDUCED)
                retval = faceYReducedIndices(i, j, k - 1);
            else if (!z_oob(i, j, k) && faceZLabels(i, j, k) == MaterialLabels::REDUCED)
                retval = faceZReducedIndices(i, j, k);
            else if (!z_oob(i, j - 1, k) && faceZLabels(i, j - 1, k) == MaterialLabels::REDUCED)
                retval = faceZReducedIndices(i, j - 1, k);
            break;
        default:
            assert(type == SamplingType::EDGEXY || type == SamplingType::EDGEXZ || type == SamplingType::EDGEYZ);
            break;
        }
    */
    return retval;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// called in construct____ActiveIndices()
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

exint
HDK_PolyStokes::Solver::serialAssignFieldIndices(
    SIM_RawIndexField& indexField,
    SIM_RawIndexField& labels,
    exint startIndex
)
{
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt* boss = UTgetInterrupt();

    exint idx = 0;

    UT_VoxelArrayIteratorI vit(indexField.fieldNC());
    //vit.setCompressOnExit(true);
    UT_VoxelTileIteratorI vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        vitt.setTile(vit);
        for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
        {
            //UT_Vector3I entry(vitt.x(), vitt.y(), vitt.z());
            if (isActive(labels.fieldNC()->getValue(vitt.x(), vitt.y(), vitt.z())))
                vitt.setValue(idx++);
        }
    }

    return idx;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// utils
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename Grid>
void
HDK_PolyStokes::Solver::uncompressTiles(Grid& grid, const UT_Array<bool>& isTileOccupiedList)
{
    UT_Interrupt* boss = UTgetInterrupt();

    UTparallelFor(UT_BlockedRange<exint>(0, isTileOccupiedList.size()), [&](const UT_BlockedRange<exint>& range)
        {
            if (boss->opInterrupt())
                return;

            for (exint tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
            {
                if (isTileOccupiedList[tileNumber])
                    grid.field()->getLinearTile(tileNumber)->uncompress();
            }
        });
}
