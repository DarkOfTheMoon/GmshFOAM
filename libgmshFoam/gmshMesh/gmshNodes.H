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
    gmshNodes

Description
    The class holds the $Nodes section of a Gmsh .msh file.

SourceFiles
    gmshNodes.C

\*---------------------------------------------------------------------------*/

#ifndef gmshNodes_H
#define gmshNodes_H

#include "IFstream.H"
#ifdef cygwin
#include "Map.hh"
#else
#include "Map.H"
#endif
#include "pointField.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// The $Nodes section data structure
class gmshNodes
{
    // Private data

    // Verbosity
    const label verbosity_;

    // Point coordinates
    pointField points_;
    // Gmsh node numbers -> OpenFOAM point numbers
    Map<label> mshToFoam_;

    // for Gmsh direct interface
    label pointI_;

public:

    gmshNodes(const label verbosity): verbosity_(verbosity) {}
    const pointField& points() const { return points_; }
    pointField& points() { return points_; }
    const Map<label>& mshToFoam() const { return mshToFoam_; }
    Map<label>& mshToFoam() { return mshToFoam_; }
    void readNodes(IFstream& inFile, const scalar formatVersion,
    const word& execName);

    // for Gmsh direct interface
    void setNVerts(const label nVerts);
    void insertNode(const label mshLabel, const scalar xVal, const scalar yVal,
    const scalar zVal);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
