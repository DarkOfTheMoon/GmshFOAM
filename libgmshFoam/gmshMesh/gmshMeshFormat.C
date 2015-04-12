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
    Reads the $MeshFormat section of a .msh file.

\*---------------------------------------------------------------------------*/

#include "gmshMeshFormat.H"
#include "gmshMessageStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gmshMeshFormat::readMeshFormat(IFstream& inFile, const word& execName)
{
    label fileType, dataSize;

    inFile >> formatVersion_ >> fileType >> dataSize;

    if(formatVersion_ != 2.0)
    {
        FatalErrorIn(execName)
            << "Unsupported $MeshFormat version " << formatVersion_
                << exit(FatalError);
    }

    gInfo(verbosity_ >= 3) << "File format version: " << formatVersion_
        << endl;

    word tag(inFile);
    if(tag != "$EndMeshFormat" && tag != "$EndMeshFormat\r")
    {
        FatalErrorIn(execName)
            << "Did not find $EndMeshFormat tag on line "
                << inFile.lineNumber() << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
