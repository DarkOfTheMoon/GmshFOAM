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
    gmshMesh

Description
    The class holds all the data structure of a Gmsh .msh file.

SourceFiles
    gmshMesh.C

\*---------------------------------------------------------------------------*/

#ifndef gmshMesh_H
#define gmshMesh_H

#include "gmshMeshFormat.H"
#include "gmshPhysicalNames.H"
#include "gmshNodes.H"
#include "gmshElements.H"

namespace Foam
{

class gmshMesh
:
    public gmshMeshFormat,
    public gmshPhysicalNames,
    public gmshNodes,
    public gmshElements
{
    // Private data

    // Verbosity
    const label verbosity_;

    // Private functions
    bool correctOrientation(const pointField& points, const cellShape& shape)
      const;
    void readSections(IFstream& inFile, const word& execName);

public:

    // Constructor
    gmshMesh(const label verbosity);
    gmshMesh(const fileName& mshName, const word& execName,
    const label verbosity);

    // Member Function
    void doAutoInvert();
    void doUnusedPointRemoval();
    void patchAndZoneNames(wordList& boundaryPatchNames,
    wordList& boundaryPatchTypes, wordList& cellZoneNames);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
