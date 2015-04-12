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

Class
    polyMeshConversion

Description
    The class holds the information and functions necessary to convert
    a gmsh data structure to polyMesh and to renumber the mesh in such
    a way that the band of the matrix is reduced. The algorithm uses a
    simple search through the neighbour list.

SourceFiles
    polyMeshConversion.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshConversion_H
#define polyMeshConversion_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class polyMeshConversion Declaration
\*---------------------------------------------------------------------------*/

class polyMeshConversion
:
    public polyMesh
{
    // Private Member Functions

    //- Disallow default bitwise copy construct
    polyMeshConversion(const polyMeshConversion&);

    //- Disallow default bitwise assignment
        void operator=(const polyMeshConversion&);

    //- Return pointer to renumbered mesh created on freestore
    polyMeshConversion* renumberedMesh(const label verbosity) const;

    label findFace(const primitivePatch&, const labelList&, const Map<label>&)
        const;

    label findInternalFace(const labelList&) const;

public:

    // Constructors

    //- Construct from points, faces and cells
    polyMeshConversion
    (
        const IOobject& io, const pointField& points, const faceList& faces,
        const cellList& cells
    );

    //- Construct from cellShapeList
    polyMeshConversion
    (
        const IOobject& io, const pointField& points,
        const cellShapeList& shapes, const faceListList& boundaryFaces,
        const wordList& boundaryPatchNames, const wordList& boundaryPatchTypes,
        const word& defaultBoundaryPatchType,
        const wordList& boundaryPatchPhysicalTypes
    );

    // Member Functions

    void repatch
    (
        List<DynamicList<label> >& zoneFaces,
        const List<DynamicList<face> >& patchFaces,
        const label verbosity
    );
    void printPatchZoneToStr(const label verbosity) const;
    void removeEmptyPatches(const label verbosity);
    void constructZones
    (
        const List<DynamicList<label> >& zoneCells,
        const List<DynamicList<label> >& zoneFaces,
        const wordList& cellZoneNames,
        const wordList& boundaryPatchNames,
        const DynamicList<label>& patchToRegion,
        const labelList& physicalNumbers
    );
    polyMeshConversion* bandCompressedMesh(const label verbosity);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
