#include "GridClassifier.h"
// houdini
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <PRM/PRM_Include.h>
#include <SIM/SIM_Object.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_DopDescription.h>
#include <UT/UT_PerfMonAutoEvent.h>
#include <SIM/SIM_FieldUtils.h>
#include <SIM/SIM_VolumetricConnectedComponentBuilder.h>
#include <UT/UT_ThreadedAlgorithm.h>
#include <UT/UT_ParallelUtil.h>
#include <tbb/tbb.h>

void
GridClassifier::classifyCells()
{
    // in: centerLabels = MaterialLabels::UNASSIGNED
    // out: centerLabels = MaterialLabels::{SOLID, GENERICFLUID, UNSOLVED}
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;
    using SIM::FieldUtils::cellToFaceMap;

    UT_Interrupt* boss = UTgetInterrupt();

    UTparallelForEachNumber(mySurfaceFieldData.field()->numTiles(), [&](const UT_BlockedRange<int>& range)
        {
            UT_VoxelArrayIteratorF vit;
            vit.setConstArray(mySurfaceFieldData.field());
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
                        if (getFieldValue(*liquidWeights.center(), cell) > 0.)
                            isInSolve = true;

                        // if any adjacent faces are inSolve, then cell should also be inSolve
                        if (!isInSolve)
                            for (int axis : {0, 1, 2})
                                for (int direction : {0, 1})
                                {
                                    UT_Vector3I face = cellToFaceMap(cell, axis, direction);
                                    if (getFieldValue(*liquidWeights.face(axis), face) > 0.)
                                        isInSolve = true;
                                }

                        if (getFieldValue(*fluidWeights.center(), cell) == 0.)
                            isInFluid = false;

                        if (isInSolve)
                        {
                            if (!isInFluid)
                                setFieldValue(*labels.center(), cell, (int)MaterialLabels::SOLID);
                            else
                                setFieldValue(*labels.center(), cell, (int)MaterialLabels::GENERICFLUID);
                        }
                        else
                        {
                            setFieldValue(*labels.center(), cell, (int)MaterialLabels::UNSOLVED);
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
GridClassifier::constructReducedRegions()
{
    // in: centerLabels = MaterialLabels::{SOLID, GENERICFLUID, UNSOLVED}
    // out: centerLabels = MaterialLabels::{SOLID, ACTIVEFLUID, REDUCED, UNSOLVED}

    constructAirBoundaryLayer();
    constructSolidBoundaryLayer();
    if (myDoTile)
        constructTiles();
    overwriteIndices(centerLabels, MaterialLabels::GENERICFLUID, MaterialLabels::REDUCED);
}
