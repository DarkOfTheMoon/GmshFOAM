/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C)
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
    Converts foam polyMesh to Gmsh mesh.

\*---------------------------------------------------------------------------*/

#include "gmshFoamConfig.H"

#include "argList.H"

#include "gmshMessageStream.H"
#include "polyMeshToGmsh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:
int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("verbosity",
    "verbosity (0-5; defaults to 3)");

#   include "setRootCase.H"

    gInfo << endl;

#if WITH_NOREAD_TIME_CTOR
    // with this Time object mesh can be written without system/ subdirectory
    Time runTime(args.rootPath(), args.caseName());
#else
#   include "createTime.H"
#endif

    label verbosity = 3; // along with the gmsh default setting of 3
    if(args.options().found("verbosity"))
    {
        verbosity = readLabel(IStringStream(args.options()["verbosity"])());
    }

    gInfo(verbosity >= 3) << "Create mesh\n" << endl;
    polyMeshToGmsh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        verbosity
    );

    mesh.writeGmshMesh();

    gInfo(verbosity >= 1) << endl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

