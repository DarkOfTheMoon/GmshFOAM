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
    gmshPhysicalNames

Description
    The class holds the $PhysicalNames section of a Gmsh .msh file.

SourceFiles
    gmshPhysicalNames.C

\*---------------------------------------------------------------------------*/

#ifndef gmshPhysicalNames_H
#define gmshPhysicalNames_H

#include "IFstream.H"
#include "labelList.H"
#include "stringList.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// The $PhysicalNames section data structure
class gmshPhysicalNames
{
    // Private data

    // Verbosity
    const label verbosity_;

    // Physical region numbers corresponding to string labels
    labelList physicalNumbers_;
    // The string labels
    stringList physicalNames_;

    label nameI_;

public:

    gmshPhysicalNames(const label verbosity)
        : verbosity_(verbosity), physicalNumbers_(0), physicalNames_(0) {}
    const labelList& physicalNumbers() const { return physicalNumbers_; }
    const stringList& physicalNames() const { return physicalNames_; }

    void readPhysicalNames(IFstream& inFile, const word& execName);

    // for Gmsh direct interface
    void setNPhysicalNames(const label nPhysicalNames);
    void insertPhysicalNames(const label physicalNumber,
    const char* physicalName);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
