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
    gmshMeshFormat

Description
    The class holds the $MeshFormat section of a Gmsh .msh file.

SourceFiles
    gmshMeshFormat.C

\*---------------------------------------------------------------------------*/

#ifndef gmshMeshFormat_H
#define gmshMeshFormat_H

#include "IFstream.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class gmshMeshFormat
{
    // Private data

    // Verbosity
    const label verbosity_;

    // The detected version of the .msh file format
    scalar formatVersion_;

public:
    gmshMeshFormat(const label verbosity): verbosity_(verbosity) {}
    scalar& formatVersion() { return formatVersion_; }
    void readMeshFormat(IFstream& inFile, const word& execName);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
