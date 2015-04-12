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
    Reads the $Nodes section of a .msh file.

\*---------------------------------------------------------------------------*/

#include "gmshNodes.H"
#include "gmshMessageStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gmshNodes::setNVerts(const label nVerts)
{
    points_.setSize(nVerts);
    pointI_ = 0;
}

void gmshNodes::insertNode(const label mshLabel, const scalar xVal,
const scalar yVal, const scalar zVal)
{
    point& pt = points_[pointI_];

    pt.x() = xVal;
    pt.y() = yVal;
    pt.z() = zVal;

    mshToFoam_.insert(mshLabel, pointI_);

    pointI_++;
}

void gmshNodes::readNodes(IFstream& inFile, const scalar formatVersion,
const word& execName)
{
    label nVerts;
    inFile >> nVerts;

    gInfo(verbosity_ >= 3) << endl << "Read nVerts: " << nVerts << endl;

    setNVerts(nVerts);

    forAll(points(), pointI)
    {
        label mshLabel;
        scalar xVal, yVal, zVal;

        inFile >> mshLabel >> xVal >> yVal >> zVal;

        insertNode(mshLabel, xVal, yVal, zVal);
    }

    word tag(inFile);
    if ((formatVersion == 1.0 && tag != "$ENDNOD" && tag != "$ENDNOD\r")
    || (formatVersion == 2.0 && tag != "$EndNodes" && tag != "$EndNodes\r"))
    {
        FatalErrorIn(execName)
            << "Did not find $ENDNOD nor $EndNodes tags on line "
                << inFile.lineNumber() << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
