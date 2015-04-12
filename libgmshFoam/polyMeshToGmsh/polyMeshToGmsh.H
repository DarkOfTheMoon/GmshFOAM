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
    Converts polyMesh to Gmsh mesh structure.

SourceFiles
    polyMeshToGmsh.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshToGmsh_H
#define polyMeshToGmsh_H

#include "polyMesh.H"
#include "cellModeller.H"
#include "Time.H"

namespace Foam
{

class polyMeshToGmsh
:
    public polyMesh
{
    // Private data

    // Verbosity
    const label verbosity_;

    // number of zone faces
    label nZoneFaces_, faceZoneI_, zoneFaceI_;

    // current physical name
    label physicalNameI_;

    // current vertex
    label vertexI_;

    // current element
    label elementI_;

    // processor No
    label myProcNo_;

    // processor No converted to string
    string myProcNoStr_;

    // number of total physical groups
    static label totalPhysicals_;

public:

    // Constructor
    polyMeshToGmsh(const IOobject& io, const label verbosity);

    // Destructor
    ~polyMeshToGmsh() { totalPhysicals_ += physicalNameI_; }

    label getNNames();
    label getNVerts() { return nPoints(); }
    label getNElems();
    void getPhysicalName(int& num, char name[]);
    void getVertex(int& num, double& x, double& y, double& z);
    void getElementAttributes(int& num, int& type, int& physical,
    int& elementary, int& partition) const;
    void getElementVerticesIndices(int indices[]);
    void writeGmshMesh();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
