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
    gmshViewSurfaceMeshMotion

Description
    Reads an OpenFOAM surfaceMesh motion and converts to a Gmsh View.

SourceFiles
    gmshViewSurfaceMeshMotion.C

\*---------------------------------------------------------------------------*/

#ifndef gmshViewSurfaceMeshMotion_H
#define gmshViewSurfaceMeshMotion_H

#include "gmshViewSurfaceMesh.H"

namespace Foam
{

template <class T1>
class gmshViewSurfaceMeshMotion
: public gmshViewSurfaceMesh<T1, vector>
{
    // Reference to the mesh
    const fvMesh& mesh_;

    // Reference to the views
    const gmshViews& views_;

    // Reference to the field type
    const label& fieldType_;

    // Reference to the runTime object
    Time& runTime_;

    // Verbosity
    const label verbosity_;

    void writeFaceMotions(char *values, const DynamicList<label>& shapeFaces,
    pointField& meshMotion, const label timeI) const;
    void writePatchMotions(char *values,
    const List<DynamicList<label> >& shapeFaces, const pointField& meshMotion,
    const label timeI, const label faceOffset) const;

public:
    // Constructor
    gmshViewSurfaceMeshMotion(const gmshViews& views,
    const word& fieldName, Time& runTime, const fvMesh& mesh,
    const label verbosity);

    // Member functions
    virtual void getViewData(gmshViewBase::gmshViewBuffer& vB) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
