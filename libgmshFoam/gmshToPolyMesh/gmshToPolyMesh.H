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

SourceFiles
    gmshToPolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef gmshToPolyMesh_H
#define gmshToPolyMesh_H

#include "Time.H"
#include "gmshMesh.H"
#include "polyMeshConversion.H"

namespace Foam
{

struct gmshToPolyMeshOptions
{
    label nOptions_;
    bool autoInvert_;
    bool checkMesh_;
    bool renumberMesh_;
    bool removeUnusedPoints_;
    label verbosity_;
};

polyMeshConversion* gmshToPolyMesh(gmshMesh& gmsh, const Time& runTime,
const gmshToPolyMeshOptions& opt);

bool writePolyMeshWithSets(const polyMeshConversion& mesh, Time& runTime,
const gmshToPolyMeshOptions& opt);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
