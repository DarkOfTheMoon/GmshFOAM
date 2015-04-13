/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Create polyMesh from cell and patch shapes.

    Suppressed a warning about undefined faces and removed unnecessary
    codes from the original polyMeshFromShapeMesh.C.

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"
#include "Time.H"
#include "primitiveMesh.H"
#include "DynamicList.H"

#include "gmshMessageStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Construct from cell shapes
polyMesh::polyMesh
(
    const IOobject& io,
    const pointField& points,
    const cellShapeList& cellsAsShapes,
    const faceListList& boundaryFaces,
    const wordList& boundaryPatchNames,
    const wordList& boundaryPatchTypes,
    const word& defaultBoundaryPatchType,
    const wordList& boundaryPatchPhysicalTypes
)
:
    objectRegistry(io),
    primitiveMesh(),
    points_
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        points
    ),
    faces_
    (
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    allOwner_
    (
        IOobject
        (
            "owner",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    allNeighbour_
    (
        IOobject
        (
            "neighbour",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    boundary_
    (
        IOobject
        (
            "boundary",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        boundaryFaces.size() + 1    // add room for a default patch
    ),
    pointZones_
    (
        IOobject
        (
            "pointZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    faceZones_
    (
        IOobject
        (
            "faceZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    cellZones_
    (
        IOobject
        (
            "cellZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    parallelDataPtr_(NULL),
    moving_(false),
    curMotionTimeIndex_(time().timeIndex()),
    oldPointsPtr_(NULL)
{
    if (debug)
    {
        Info <<"Constructing polyMesh from cell and boundary shapes." << endl;
    }

    // Remove all of the old mesh files if they exist
    removeFiles(instance());

    // Calculate the faces of all cells
    // Initialise maximum possible numer of mesh faces to 0
    label maxFaces = 0;

    // Set up a list of face shapes for each cell
    faceListList cellsFaceShapes(cellsAsShapes.size());
    cellList cells(cellsAsShapes.size());

    forAll(cellsFaceShapes, cellI)
    {
        cellsFaceShapes[cellI] = cellsAsShapes[cellI].faces();

        cells[cellI].setSize(cellsFaceShapes[cellI].size());

        // Initialise cells to -1 to flag undefined faces
        static_cast<labelList&>(cells[cellI]) = -1;

        // Count maximum possible numer of mesh faces
        maxFaces += cellsFaceShapes[cellI].size();
    }

    // Set size of faces array to maximum possible number of mesh faces
    faces_.setSize(maxFaces);

    // Initialise number of faces to 0
    label nFaces = 0;

    // set reference to point-cell addressing
    labelListList PointCells = cellShapePointCells(cellsAsShapes);

    bool found = false;

    forAll(cells, cellI)
    {
        // Note:
        // Insertion cannot be done in one go as the faces need to be
        // added into the list in the increasing order of neighbour
        // cells.  Therefore, all neighbours will be detected first
        // and then added in the correct order.  

        const faceList& curFaces = cellsFaceShapes[cellI];

        // Record the neighbour cell
        labelList neiCells(curFaces.size(), -1);

        // Record the face of neighbour cell
        labelList faceOfNeiCell(curFaces.size(), -1);

        label nNeighbours = 0;

        // For all faces ...
        forAll(curFaces, faceI)
        {
            // Skip faces that have already been matched
            if (cells[cellI][faceI] >= 0) continue;

            found = false;

            const face& curFace = curFaces[faceI];

            // Get the list of labels
            const labelList& curPoints = curFace;

            // For all points
            forAll(curPoints, pointI)
            {
                // dGget the list of cells sharing this point
                const labelList& curNeighbours =
                    PointCells[curPoints[pointI]];

                // For all neighbours
                forAll(curNeighbours, neiI)
                {
                    label curNei = curNeighbours[neiI];

                    // Reject neighbours with the lower label
                    if (curNei > cellI)
                    {
                        // Get the list of search faces
                        const faceList& searchFaces = cellsFaceShapes[curNei];

                        forAll(searchFaces, neiFaceI)
                        {
                            if (searchFaces[neiFaceI] == curFace)
                            {
                                // Match!!
                                found = true;

                                // Record the neighbour cell and face
                                neiCells[faceI] = curNei;
                                faceOfNeiCell[faceI] = neiFaceI;
                                nNeighbours++;

                                break;
                            }
                        }
                        if (found) break;
                    }
                    if (found) break;
                }
                if (found) break;
            } // End of current points
        }  // End of current faces

        // Add the faces in the increasing order of neighbours
        for (label neiSearch = 0; neiSearch < nNeighbours; neiSearch++)
        {
            // Find the lowest neighbour which is still valid
            label nextNei = -1;
            label minNei = cells.size();

            forAll (neiCells, ncI)
            {
                if (neiCells[ncI] > -1 && neiCells[ncI] < minNei)
                {
                    nextNei = ncI;
                    minNei = neiCells[ncI];
                }
            }

            if (nextNei > -1)
            {
                // Add the face to the list of faces
                faces_[nFaces] = curFaces[nextNei];

                // Set cell-face and cell-neighbour-face to current face label
                cells[cellI][nextNei] = nFaces;
                cells[neiCells[nextNei]][faceOfNeiCell[nextNei]] = nFaces;

                // Stop the neighbour from being used again
                neiCells[nextNei] = -1;

                // Increment number of faces counter
                nFaces++;
            }
            else
            {
                FatalErrorIn
                (
                    "polyMesh::polyMesh\n"
                    "(\n"
                    "    const IOobject& io,\n"
                    "    const pointField& points,\n"
                    "    const cellShapeList& cellsAsShapes,\n"
                    "    const faceListList& boundaryFaces,\n"
                    "    const wordList& boundaryPatchTypes,\n"
                    "    const wordList& boundaryPatchNames,\n"
                    "    const word& defaultBoundaryPatchType\n"
                    ")"
                )   << "Error in internal face insertion"
                    << abort(FatalError);
            }
        }
    }

    // Do boundary faces

    labelList patchSizes(boundaryFaces.size(), 0);
    labelList patchStarts(boundaryFaces.size(), nFaces);

    // Grab "non-existing" faces and put them into a default patch

    label defaultPatchStart = nFaces;

    forAll(cells, cellI)
    {
        labelList& curCellFaces = cells[cellI];

        forAll(curCellFaces, faceI)
        {
            if (curCellFaces[faceI] == -1) // "non-existent" face
            {
                curCellFaces[faceI] = nFaces;
                faces_[nFaces] = cellsFaceShapes[cellI][faceI];

                nFaces++;
            }
        }
    }

    // Reset the size of the face list
    faces_.setSize(nFaces);

    // Warning: Patches can only be added once the face list is
    // completed, as they hold a subList of the face list
    forAll (boundaryFaces, patchI)
    {
        // add the patch to the list
        boundary_.hook
        (
            polyPatch::New
            (
                boundaryPatchTypes[patchI],
                boundaryPatchNames[patchI],
                patchSizes[patchI],
                patchStarts[patchI],
                patchI,
                boundary_
            )
        );

        if
        (
            boundaryPatchPhysicalTypes.size()
         && boundaryPatchPhysicalTypes[patchI].size()
        )
        {
            boundary_[patchI].physicalType() =
                boundaryPatchPhysicalTypes[patchI];
        }
    }

    label nAllPatches = boundaryFaces.size();

    if (nFaces > defaultPatchStart)
    {
        nAllPatches++;

        boundary_.hook
        (
            polyPatch::New
            (
                defaultBoundaryPatchType,
                "defaultFaces",
                nFaces - defaultPatchStart,
                defaultPatchStart,
                boundary_.size() - 1,
                boundary_
            )
        );
    }

    // Reset the size of the boundary
    boundary_.setSize(nAllPatches);

    // Set the primitive mesh
    initMesh(cells);

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
