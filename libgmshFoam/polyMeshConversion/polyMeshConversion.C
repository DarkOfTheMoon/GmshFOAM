/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
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
    The class holds the information and functions necessary to convert
    a gmsh data structure to polyMesh and to renumber the mesh in such
    a way that the band of the matrix is reduced. The algorithm is
    cited from the fvMeshBandCompression class in the renumberMesh
    utility.

    2) Added cellZone/faceZone handlings. 2007/3/11 Takuya OSHIMA

\*---------------------------------------------------------------------------*/

#include "bandCompression.H"
#include "SLList.H"
#include "Time.H"
#include "SortableList.H"
#include "repatchPolyTopoChanger.H"

#include "gmshMessageStream.H"
#include "polyMeshConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

polyMeshConversion::polyMeshConversion
(
    const IOobject& io,
    const Xfer<pointField>& points,
    const Xfer<faceList>& faces,
    const Xfer<cellList>& cells
)
    :
    polyMesh(io, points, faces, cells)
{}

polyMeshConversion::polyMeshConversion
(
    const IOobject& io,
    const Xfer<pointField>& points,
    const cellShapeList& shapes,
    const faceListList& boundaryFaces,
    const wordList& boundaryPatchNames,
    const wordList& boundaryPatchTypes,
    const word& defaultBoundaryPatchName,
    const word& defaultBoundaryPatchType,
    const wordList& boundaryPatchPhysicalTypes
)
:
    polyMesh(io, points, shapes, boundaryFaces, boundaryPatchNames,
    boundaryPatchTypes,defaultBoundaryPatchName, defaultBoundaryPatchType, boundaryPatchPhysicalTypes)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

polyMeshConversion* polyMeshConversion::renumberedMesh(const label verbosity)
    const
{

    //- New cell order for the compressed band
    // This gives the order in which the cells in the original mesh
    // are to be addressed to create a mesh with a compressed band
    // The order list gives the old cell label for every new cell
    const labelList order(bandCompression(cellCells()));

    // Renumber the cell list

    const cellList& oldCells = cells();
    const labelList& oldOwner = faceOwner();
    const labelList& oldNeighbour = faceNeighbour();

    cellList newCells(oldCells.size());

    // The reverse order list gives the new cell label for every old cell
    labelList reverseOrder(order.size());

    forAll (order, cellI)
    {
        newCells[cellI] = oldCells[order[cellI]];

        reverseOrder[order[cellI]] = cellI;
    }

    // Renumber the faces.
    // Reverse face order gives the new face number for every old face
    labelList reverseFaceOrder(nFaces(), 0);

    // Mark the internal faces with -2 so that they are inserted first
    forAll (newCells, cellI)
    {
        const labelList& curFaces = newCells[cellI];

        forAll (curFaces, faceI)
        {
            reverseFaceOrder[curFaces[faceI]]--;
        }
    }

    // Order internal faces
    label nMarkedFaces = 0;

    forAll (newCells, cellI)
    {
        // Note:
        // Insertion cannot be done in one go as the faces need to be
        // added into the list in the increasing order of neighbour
        // cells.  Therefore, all neighbours will be detected first
        // and then added in the correct order.  

        const labelList& curFaces = newCells[cellI];

        // Record the neighbour cell
        labelList neiCells(curFaces.size(), -1);

        label nNeighbours = 0;

        forAll (curFaces, faceI)
        {
            if (reverseFaceOrder[curFaces[faceI]] == -2)
            {
                // Face is internal and gets reordered
                if (cellI == reverseOrder[oldOwner[curFaces[faceI]]])
                {
                    neiCells[faceI] =
                        reverseOrder[oldNeighbour[curFaces[faceI]]];
                }
                else if (cellI == reverseOrder[oldNeighbour[curFaces[faceI]]])
                {
                    neiCells[faceI] =
                        reverseOrder[oldOwner[curFaces[faceI]]];
                }
                else
                {
                    gWarning(verbosity >= 1)
                        << "Screwed up in renumbering mesh!!!" << endl;
                }

                nNeighbours++;
            }
        }

        // Add the faces in the increasing order of neighbours
        for (label neiSearch = 0; neiSearch < nNeighbours; neiSearch++)
        {
            // Find the lowest neighbour which is still valid
            label nextNei = -1;
            label minNei = oldCells.size();

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
                // Face is internal and gets reordered
                reverseFaceOrder[curFaces[nextNei]] = nMarkedFaces;

                // Stop the neighbour from being used again
                neiCells[nextNei] = -1;

                nMarkedFaces++;
            }
            else
            {
                gWarning(verbosity >= 1)
                    << "Error in internal face insertion; " << nl
                    << "    falling back to the mesh without matrix bandwidth"
                    << " compression" << endl;
                return static_cast<polyMeshConversion*>(NULL);
            }
        }
    }
    // Insert the boundary faces into reordering list
    forAll (reverseFaceOrder, faceI)
    {
        if (reverseFaceOrder[faceI] < 0)
        {
            reverseFaceOrder[faceI] = nMarkedFaces;

            nMarkedFaces++;
        }
    }

    // Face order gives the old face label for every new face
    labelList faceOrder(reverseFaceOrder.size());

    forAll (faceOrder, faceI)
    {
        faceOrder[reverseFaceOrder[faceI]] = faceI;
    }

    // Renumber the cells
    forAll (newCells, cellI)
    {
        labelList oldF = newCells[cellI];

        labelList& newF = newCells[cellI];

        forAll (newF, fI)
        {
            newF[fI] = reverseFaceOrder[oldF[fI]];
        }
    }

    faceList newFaces(this->faces().size());

    const faceList& oldFaces = this->faces();

    forAll (newFaces, faceI)
    {
        newFaces[faceI] = oldFaces[faceOrder[faceI]];
    }

    // Turn the face that need to be turned

    // Only loop through internal faces
    forAll (oldNeighbour, faceI)
    {
        if
        (
            reverseOrder[oldNeighbour[faceOrder[faceI]]]
          < reverseOrder[oldOwner[faceOrder[faceI]]]
        )
        {
            newFaces[faceI] = newFaces[faceI].reverseFace();
        }
    }

    // Make a new mesh
    polyMeshConversion* newMeshPtr = new polyMeshConversion
    (
        IOobject
        (
            polyMesh::defaultRegion,
            time().constant(),
            time()
        ),
        Xfer<pointField>(points()),
        Xfer<faceList>(newFaces),
        Xfer<cellList>(newCells)
    );

    polyMeshConversion& newMesh = *newMeshPtr;

    // Add boundaries

    const polyPatchList& oldPatches = boundaryMesh();

    List<polyPatch*> newPatches(oldPatches.size());

    forAll (newPatches, patchI)
    {
        newPatches[patchI] =
            oldPatches[patchI].clone(newMesh.boundaryMesh()).ptr();
    }

    newMesh.addPatches(newPatches);

    // Add zones

    List<cellZone*> cz(cellZones().size());
    forAll(cellZones(), zoneI)
    {
        SortableList<label> zoneCells(cellZones()[zoneI].size());

        forAll(cellZones()[zoneI], cellI)
        {
            zoneCells[cellI] = reverseOrder[cellZones()[zoneI][cellI]];
        }
	zoneCells.sort();

        cz[zoneI] = new cellZone(cellZones().names()[zoneI], zoneCells, zoneI,
				 newMesh.cellZones());
    }

    List<faceZone*> fz(faceZones().size());
    forAll(faceZones(), zoneI)
    {
        SortableList<label> zoneFaces(faceZones()[zoneI].size());

        forAll(faceZones()[zoneI], faceI)
        {
            zoneFaces[faceI] = reverseFaceOrder[faceZones()[zoneI][faceI]];
        }
	zoneFaces.sort();

        fz[zoneI] = new faceZone(faceZones().names()[zoneI], zoneFaces,
                boolList(zoneFaces.size(), true), zoneI, newMesh.faceZones());
    }

    if (cz.size() > 0 || fz.size() > 0)
    {
        newMesh.addZones(List<pointZone*>(0), fz, cz);
    }

    return newMeshPtr;
}

// Find face in pp which uses all vertices in meshF (in mesh point labels)
label polyMeshConversion::findFace(const primitivePatch& pp,
const labelList& meshF, const Map<label> &meshPointMap) const
{
    // move creation of map to outer loop
    // const Map<label>& meshPointMap = pp.meshPointMap();

    // meshF[0] in pp labels.
    if (!meshPointMap.found(meshF[0]))
    {
        // not issue a warning because the undetermined face will successively
        // be handled by findInternalFace.
        // gWarning<< "Not using gmsh face " << meshF
        //     << " since zero vertex is not on boundary of polyMesh" << endl;
        return -1;
    }

    // Find faces using first point
    const labelList& pFaces = pp.pointFaces()[meshPointMap[meshF[0]]];

    // Go through all these faces and check if there is one which uses all of
    // meshF vertices (in any order ;-)
    forAll(pFaces, i)
    {
        label faceI = pFaces[i];

        const face& f = pp[faceI];

        // Count uses of vertices of meshF for f
        label nMatched = 0;

        forAll(f, fp)
        {
            if (findIndex(meshF, f[fp]) != -1)
            {
                nMatched++;
            }
        }

        if (nMatched == meshF.size())
        {
            return faceI;
        }
    }

    return -1;
}


// Same but find internal face. Expensive addressing.
label polyMeshConversion::findInternalFace(const labelList& meshF) const
{
    const labelList& pFaces = pointFaces()[meshF[0]];

    forAll(pFaces, i)
    {
        label faceI = pFaces[i];

        const face& f = faces()[faceI];

        // Count uses of vertices of meshF for f
        label nMatched = 0;

        forAll(f, fp)
        {
            if (findIndex(meshF, f[fp]) != -1)
            {
                nMatched++;
            }
        }

        if (nMatched == meshF.size())
        {
            return faceI;
        }        
    }
    return -1;
}

void polyMeshConversion::repatch(List<DynamicList<label> >& zoneFaces,
const List<DynamicList<face> >& patchFaces, const label verbosity)
{
    repatchPolyTopoChanger repatcher(*this);

    // Now use the patchFaces to patch up the outside faces of the mesh.

    // Get the patch for all the outside faces (= default patch added as last)
    const polyPatch& pp = boundaryMesh()[boundaryMesh().size()-1];
    const Map<label> boundaryMeshPointMap = pp.meshPointMap();

    zoneFaces.setSize(patchFaces.size());

    gInfo(verbosity >= 4) << nl
        << "Finding corresponding patch faces for the " << pp.size()
        << " undefined faces:" << endl;

    label nPatchFaces = 0, nZoneFaces = 0, nUnfoundFaces = 0;
    // Go through all the patchFaces and find corresponding face in pp.
    forAll(patchFaces, patchI)
    {
        const DynamicList<face>& pFaces = patchFaces[patchI];

        gInfo(verbosity >= 5) << "    Finding faces of "
            << boundaryMesh()[patchI].name() << endl;

        forAll(pFaces, i)
        {
            gInfo(verbosity >= 5 && (i % 1000) == 0) << "\tFinding face "
                << i << " of " << pFaces.size() << endl;

            const face& f = pFaces[i];

            // Find face in pp using all vertices of f.
            //label patchFaceI = findFace(pp, f);
            label patchFaceI = findFace(pp, f, boundaryMeshPointMap);

            if (patchFaceI != -1)
            {
                label meshFaceI = pp.start() + patchFaceI;

                repatcher.changePatchID(meshFaceI, patchI);
                nPatchFaces++;
            }
            else
            {
                // Maybe internal face? If so add to faceZone with same index
                // - might be useful.
                label meshFaceI = findInternalFace(f);

                if (meshFaceI != -1)
                {
                    zoneFaces[patchI].append(meshFaceI);
                    nZoneFaces++;
                }
                else
                {
		    nUnfoundFaces++;
		    gWarning((verbosity >= 1 && nUnfoundFaces <= 50)
                    || verbosity >= 5) << "    Could not match face "
                        << f << endl;
                }
            }
        }
    }
    if(verbosity >= 5)
    {
        gWarning(nUnfoundFaces > 0) << "        ... " << nUnfoundFaces
            << " unmatched faces in total." << endl;
    }
    else 
    {
        gWarning(verbosity >= 1 && nUnfoundFaces > 50) << "        ... "
            << nUnfoundFaces - 50 << " more unmatched faces follow ("
            << nUnfoundFaces << " in total)" << endl;
    }

    gInfo(verbosity >= 4) << "    ... defined " << nPatchFaces
        << " patch faces and " << nZoneFaces << " faceZone faces." << endl;

    Info << nl << "Performing repatching:" << endl;

    repatcher.repatch();

    Info << "    ... repatching done." << endl;

    forAll(zoneFaces, zoneI)
    {
        zoneFaces[zoneI].shrink();
    }
}

void polyMeshConversion::printPatchZoneToStr(const label verbosity) const
{
    // Print the finally determined patches
    label maxLen = 0;
    forAll(boundaryMesh(), patchI)
    {
        label len = boundaryMesh()[patchI].name().length();
        if (len > maxLen)
        {
            maxLen = len;
        }
    }

    Info << nl << "Patches:" << nl << "Patch\tSize\tName";

    for(label i = 0; i < maxLen - 4; i++)
    {
        Info << ' ';
    }
    Info << "\tBase type" << endl;

    forAll(boundaryMesh(), patchI)
    {
        label len = boundaryMesh()[patchI].name().length();
        Info << "    " << patchI << '\t'
            << boundaryMesh()[patchI].size() << '\t'
            << boundaryMesh()[patchI].name();
        for(label i = 0; i < maxLen - len; i++)
        {
            Info << ' ';
        }
        Info << '\t' << boundaryMesh()[patchI].type() << endl;
    }

    // Print the finally determined cellZones
    if (cellZones().size())
    {
        Info << nl << "CellZones:" << nl << "Zone\tSize\tName"
            << endl;

        forAll(cellZones(), zoneI)
        {
            Info << "    " << zoneI << '\t'
                << cellZones()[zoneI].size() << '\t'
                << cellZones()[zoneI].name() << endl;
        }
    }

    // Print the finally determined faceZones
    if (faceZones().size())
    {
        Info << nl << "FaceZones:" << nl << "Zone\tSize\tName"
            << endl;

        forAll(faceZones(), zoneI)
        {
            Info << "    " << zoneI << '\t'
                << faceZones()[zoneI].size() << '\t'
                << faceZones()[zoneI].name() << endl;
        }
    }
}

void polyMeshConversion::removeEmptyPatches(const label verbosity)
{
    gInfo(verbosity >= 4) << nl
        << "Removing zero-sized surface patches from polyMesh/boundary"
        << (verbosity >= 5 ? ":" : "") << endl;

    for(label patchI = boundaryMesh().size() - 1; patchI >= 0; patchI--)
    {
        if(boundaryMesh()[patchI].size() == 0)
        {
            gInfo(verbosity >= 5) << "    " << boundaryMesh()[patchI].name()
                << endl;

            for(label patchJ = patchI; patchJ < boundaryMesh().size() - 1;
		patchJ++)
	    {
                const_cast<polyBoundaryMesh &>(boundaryMesh())[patchJ]
                    = boundaryMesh()[patchJ + 1];
	    }
	    const_cast<polyBoundaryMesh &>(boundaryMesh())
                .setSize(boundaryMesh().size() - 1);
	}
    }
}

void polyMeshConversion::constructZones(
    const List<DynamicList<label> >& zoneCells,
    const List<DynamicList<label> >& zoneFaces, const wordList& cellZoneNames,
    const wordList& boundaryPatchNames,
    const DynamicList<label>& patchToRegion,
    const labelList& physicalNumbers)
{
    // Construct and add the zones. Note that cell ordering does not change
    // because of repatch() and neither does internal faces so we can
    // use the zoneCells/zoneFaces as is.
    List<cellZone*> cz;
    List<faceZone*> fz;

    label nValidCellZones = 0;
    forAll(zoneCells, zoneI)
    {
        if (zoneCells[zoneI].size() > 0)
        {
            nValidCellZones++;
        }
    }

    if (nValidCellZones > 0)
    {
        cz.setSize(nValidCellZones);

        nValidCellZones = 0;

        forAll(zoneCells, zoneI)
        {
            if (zoneCells[zoneI].size() > 0)
            {
                cz[nValidCellZones] = new cellZone
                    (
                        cellZoneNames[zoneI],
                        zoneCells[zoneI],
                        nValidCellZones,
                        cellZones()
                    );
                nValidCellZones++;
            }
        }
    }

    label nValidFaceZones = 0;
    forAll(zoneFaces, zoneI)
    {
        if (zoneFaces[zoneI].size() > 0)
        {
            nValidFaceZones++;
        }
    }
    if (nValidFaceZones > 0)
    {
        fz.setSize(nValidFaceZones);

        nValidFaceZones = 0;

        forAll(zoneFaces, zoneI)
        {
            if (zoneFaces[zoneI].size() > 0)
            {
	        word zoneName;

                // Determine whether the name for the faceZone is defined,
                // and use the name if defined
                forAll(physicalNumbers, numberI)
                {
                    if(patchToRegion[zoneI] == physicalNumbers[numberI])
                    {
                        zoneName = boundaryPatchNames[zoneI];
                    }
                }
                // If undefined give an automatially generated name
                if(zoneName.length() == 0)
                {
		    zoneName = "faceZone_" + Foam::name(patchToRegion[zoneI]);
                }
                fz[nValidFaceZones] = new faceZone
                    (
                        zoneName,
                        zoneFaces[zoneI],
                        boolList(zoneFaces[zoneI].size(), true),
                        nValidFaceZones,
                        faceZones()
                    );
                nValidFaceZones++;
            }
        }
    }

    if (cz.size() > 0 || fz.size() > 0)
    {
        addZones(List<pointZone*>(0), fz, cz);
    }
}

polyMeshConversion *polyMeshConversion::bandCompressedMesh(
    const label verbosity)
{
    const unallocLabelList& oldFownerC = faceOwner();
    const unallocLabelList& oldFneighbourC = faceNeighbour();

    label band = 0;

    forAll(oldFneighbourC, faceI)
    {
        label diff = oldFneighbourC[faceI] - oldFownerC[faceI];

        if (diff > band)
        {
            band = diff;
        }
    }

    Info << nl
        << "Performing renumberMesh matrix bandwidth compression";
    Info << ":" << nl <<
        "    Band before renumbering: " << band;
    Info << endl;

    // create a new mesh
    polyMeshConversion* newMeshPtr = renumberedMesh(verbosity);
    polyMeshConversion& newMesh
        = (newMeshPtr == static_cast<polyMesh*>(NULL) ? *this : *newMeshPtr);

    // check the new bandwidth
    band = 0;
    const unallocLabelList& newFownerC = newMesh.faceOwner();
    const unallocLabelList& newFneighbourC = newMesh.faceNeighbour();

    forAll(newFneighbourC, faceI)
    {
        label diff = newFneighbourC[faceI] - newFownerC[faceI];

        if (diff > band)
        {
            band = diff;
        }
    }

    Info << "    Band after renumbering: " << band << endl;

    return newMeshPtr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
