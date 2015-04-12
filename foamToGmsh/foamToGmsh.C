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
    Converts foam output fields to Gmsh view format files.

\*---------------------------------------------------------------------------*/

#include "argList.H"

#include "gmshMessageStream.H"
#include "gmshViews.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:
int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("splitTimeStepsByMeshMotion", "");
    argList::validOptions.insert("startTime",
    "time (\"latestTime\" for the latest time)");
    argList::validOptions.insert("verbosity",
    "verbosity (0-5; defaults to 3)");

#   include "setRootCase.H"
    gInfo << endl;
#   include "createTime.H"
#   include "createMesh.H"

    gmshViews::gmshViewsOptions opt;

    opt.nOptions_ = 3;
    opt.splitTimeStepsByMeshMotion_
        = args.options().found("splitTimeStepsByMeshMotion");

    if(args.options().found("startTime"))
    {
        opt.startTimeStr_ = args.options()["startTime"];
    }
    else
    {
        opt.startTimeStr_ = "";
    }

    opt.verbosity_ = 3; // along with the gmsh default setting of 3
    if(args.options().found("verbosity"))
    {
        opt.verbosity_
            = readLabel(IStringStream(args.options()["verbosity"])());
    }

    // initialization and validity check
    bool good;
    gmshViews views(good, args.rootPath(), args.caseName(), opt);
    if(!good)
    {
        FatalErrorIn(args.executable()) << "No valid field found."
            << exit(FatalError);
    }

    // set the folder where the view files go
    const fileName gmshViewsFolder = args.rootPath()/args.caseName()/"Gmsh";
    if(isDir(gmshViewsFolder))
    {
        gInfo(opt.verbosity_ >= 3) << endl << "Removing old Gmsh files in "
            << gmshViewsFolder << endl;
        rmDir(gmshViewsFolder);
    }
    mkDir(gmshViewsFolder);

    // do conversion
    views.convert(gmshViewsFolder);

    gInfo(opt.verbosity_ >= 1) << endl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

