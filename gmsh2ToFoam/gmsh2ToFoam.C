/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original authors
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
    Reads .msh file as written by Gmsh.

    Needs surface elements on mesh to be present and aligned with outside faces
    of the mesh. I.e. if the mesh is hexes, the outside faces need to be quads

    Note: There is something seriously wrong with the ordering written in the
    .msh file. Normal operation is to use the ordering as described
    in the manual. Use the -autoInvert to invert based on the geometry.
    Not very well tested.

    Note: The code now uses the element (cell,face) physical region id number
    to create cell zones and faces zones (similar to
    fluentMeshWithInternalFaces).

    A use of the cell zone information, is for field initialization with the
    "setFields" utility. see the classes:  topoSetSource, zoneToCell.  

    2)  The addition of the flag "-verbose" to control the
    amount of screen output

    Modifications made as gmsh2ToFoam are written in the HISTORY file.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "argList.H"

#include "gmshToPolyMesh.H"
#include "gmshMessageStream.H"

extern "C" {
#include <stdlib.h>
};

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void runGmsh(const fileName& geoName, const argList& args,
const label verbosity)
{
    // Check if the file is present
    IFstream *geoIn = new IFstream(geoName);
    if(geoIn->bad())
    {
        delete geoIn;
        FatalErrorIn(args.executable())
            << "The .geo file " << geoName << " not present or readable."
                << exit(FatalError);
    }
    delete geoIn;

    // Sanity check for the .geo file name
    forAll(geoName, charI)
    {
        const char c = geoName[charI];
        if(!isalnum(c) && c != '.' && c != '-' && c != '_' && c != '/')
        {
            FatalErrorIn(args.executable())
                << "The .geo file name " << geoName << " contains" << nl
                    << "    an invalid character '" << c << "'."
                    << exit(FatalError);
        }
    }

    const string cmdLine = "gmsh -3 '" + geoName + "'";

    gInfo(verbosity >= 2) << "A file name with a .geo extention is given;"
        << endl << "running Gmsh with " << cmdLine << ":" << endl;

    int status = system(cmdLine.c_str());

    if(status == -1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0)
    {
        FatalErrorIn(args.executable())
            << "Failed running gmsh. Status = "
                << status << exit(FatalError);
    }

    gInfo(verbosity >= 2) << "    ... running Gmsh done." << endl;
}

polyMeshConversion* convertedPolyMesh(const fileName& mshName,
const Time& runTime, const argList& args, const gmshToPolyMeshOptions& opt)
{
    gmshMesh gmsh(mshName, args.executable(), opt.verbosity_);

    polyMeshConversion* meshPtr = gmshToPolyMesh(gmsh, runTime, opt);

    return meshPtr;
}

// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append(".msh or .geo file");
    argList::validOptions.insert("noAutoInvert", "");
    argList::validOptions.insert("noCheckMesh", "");
    argList::validOptions.insert("noRenumberMesh", "");
    argList::validOptions.insert("noUnusedPointRemoval", "");
    argList::validOptions.insert("verbosity", "verbosity (0-5; defaults to 3)");

    #include "setRootCase.H"
    #include "createTime.H"

    gInfo << endl;

    gmshToPolyMeshOptions opt;
    opt.nOptions_ = 5;
    opt.autoInvert_ = !args.options().found("noAutoInvert");
    opt.checkMesh_ = !args.options().found("noCheckMesh");
    opt.renumberMesh_ = !args.options().found("noRenumberMesh");
    opt.removeUnusedPoints_ = !args.options().found("noUnusedPointRemoval");

    opt.verbosity_ = 3; // along with the gmsh default setting of 3
    if(args.options().found("verbosity"))
    {
        opt.verbosity_
            = readLabel(IStringStream(args.options()["verbosity"])());
    }

    fileName mshName(args.args()[3]);

    if(mshName.ext() == fileName("geo"))
    {
        runGmsh(mshName, args, opt.verbosity_);

        mshName = mshName.lessExt() + ".msh";
        gInfo(opt.verbosity_ >= 2) << endl << "Converting " << mshName << endl
            << endl;
    }

    polyMeshConversion* meshPtr
        = convertedPolyMesh(mshName, runTime, args, opt);

    writePolyMeshWithSets(*meshPtr, runTime, opt);

    delete meshPtr;

    gInfo(opt.verbosity_ >= 1) << endl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
