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
    gmshViewPointMesh

Description
    Reads an OpenFOAM pointField and converts to a Gmsh View.

SourceFiles
    gmshViewPointMesh.C

\*---------------------------------------------------------------------------*/

#ifndef gmshViewPointMesh_H
#define gmshViewPointMesh_H

#include "gmshViewBase.H"

namespace Foam
{

class gmshViews;

template <class T1, class T2>
class gmshViewPointMesh
: public gmshViewBase
{
    const static label geoOffset_;

    // Reference to the mesh
    const fvMesh& mesh_;

    void writePointValues(char *values, const GeometricField<T2,
    pointPatchField, pointMesh>& field, const label timeI) const;

protected:
    // Protected data
    const static label fieldType_;
    const static enum meshTypes meshType_;

    // Protected member functions
    virtual void writePointPoints(char *p) const;

public:
    // Constructor
    gmshViewPointMesh(const gmshViews& views, const word& fieldName,
    Time& runTime, const fvMesh& mesh, const label verbosity);

    // Member functions
    virtual const enum meshTypes meshType() const { return meshType_; }
    virtual void getPostFormat(gmshPostFormat& pF) const;
    virtual void getTimeStepValues(char *tB) const;
    virtual void getViewHeader(gmshViewHeader& vH) const;
    virtual void getViewData(gmshViewBuffer& vB) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
