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
    Converts a gmshMesh structure to a polyMesh structure.

\*---------------------------------------------------------------------------*/

#include "cellSet.H"
#include "faceSet.H"

#include "gmshMessageStream.H"
#include "gmshToPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshConversion* gmshToPolyMesh(gmshMesh& gmsh, const Time& runTime,
const gmshToPolyMeshOptions& opt)
{
    wordList boundaryPatchNames, boundaryPatchTypes, cellZoneNames;

    gmsh.patchAndZoneNames(boundaryPatchNames, boundaryPatchTypes,
    cellZoneNames);

    if(opt.autoInvert_)
    {
        gmsh.doAutoInvert();
    }

    if(opt.removeUnusedPoints_)
    {
        gmsh.doUnusedPointRemoval();
    }

    gInfo(opt.verbosity_ >= 3) << endl << "Constructing polyMesh" << endl;

    // Problem is that the orientation of the patchFaces does not have to
    // be consistent with the outwards orientation of the mesh faces. So
    // we have to construct the mesh in two stages:
    // 1. define mesh with all boundary faces in one patch
    // 2. use the read patchFaces to find the corresponding boundary face
    //    and repatch it.

    // Create correct number of patches
    // (but without any faces in it)
    faceListList boundaryFaces(gmsh.patchFaces().size());

    word defaultFacesType = polyPatch::typeName;
    wordList boundaryPatchPhysicalTypes
        (
            boundaryFaces.size(),
            polyPatch::typeName
        );

    polyMeshConversion* meshPtr = new polyMeshConversion
        (
            IOobject
            (
                polyMesh::defaultRegion,
                runTime.constant(),
                runTime
            ),
            gmsh.points(),
            gmsh.cells(),
            boundaryFaces,
            boundaryPatchNames,
            boundaryPatchTypes,
            defaultFacesType,
            boundaryPatchPhysicalTypes
        );
    polyMeshConversion& mesh = *meshPtr;

    // Storage for faceZones.
    List<DynamicList<label> > zoneFaces;

    mesh.repatch(zoneFaces, gmsh.patchFaces(), opt.verbosity_);
    mesh.removeEmptyPatches(opt.verbosity_);
    mesh.constructZones(gmsh.zoneCells(), zoneFaces, cellZoneNames,
    boundaryPatchNames, gmsh.patchToRegion(), gmsh.physicalNumbers());

    polyMeshConversion* newMeshPtr;

    if(opt.renumberMesh_)
    {
        newMeshPtr = meshPtr->bandCompressedMesh(opt.verbosity_);
        if(newMeshPtr != meshPtr) // success
        {
            delete meshPtr;
        }
    }
    else
    {
        newMeshPtr = meshPtr;
    }

    if(opt.checkMesh_)
    {
        gInfo(opt.verbosity_ >= 1) << endl
            << "Performing a simplified checkMesh test"
            << (opt.verbosity_ >= 3 ? ":" : "") << endl;
        newMeshPtr->checkMesh(opt.verbosity_ >= 4);
        gInfo(opt.verbosity_ >= 3) << "    ... checkMesh done." << endl;
    }

    newMeshPtr->printPatchZoneToStr(opt.verbosity_);

    return newMeshPtr;
}

bool writePolyMeshWithSets(const polyMeshConversion& mesh, Time& runTime,
const gmshToPolyMeshOptions& opt)
{
    gInfo(opt.verbosity_ >= 3) << endl
        << "Writing the final converted mesh and sets:" << endl;

    //Get polyMesh to write to constant
    runTime.setTime(instant(runTime.constant()), 0);

    try
    {
        mesh.write();

        if(mesh.cellZones().size() > 0)
        {
            // Write the cellSets which has the same content as cellZones
            forAll(mesh.cellZones(), zoneI)
                {
                    cellSet cset(mesh, mesh.cellZones()[zoneI].name(),
                    mesh.cellZones()[zoneI]);
                    cset.write();
                }
        }

        if(mesh.faceZones().size() > 0)
        {
            // Write the faceSets which has the same content as faceZones
            forAll(mesh.faceZones(), zoneI)
                {
                    faceSet fset(mesh, mesh.faceZones()[zoneI].name(),
                    mesh.faceZones()[zoneI]);
                    fset.write();
                }
        }
    }
    catch(error& e)
    {
        gSeriousError(opt.verbosity_ >= 1) << e.message().c_str() << endl;
        return false;
    }

    gInfo(opt.verbosity_ >= 3) << "    ... done." << endl;
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
