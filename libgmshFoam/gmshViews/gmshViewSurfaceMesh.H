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
    gmshViewSurfaceMesh

Description
    Reads an OpenFOAM surfaceField and converts to a Gmsh View.

SourceFiles
    gmshViewSurfaceMesh.C

\*---------------------------------------------------------------------------*/

#ifndef gmshViewSurfaceMesh_H
#define gmshViewSurfaceMesh_H

#include "gmshFoamConfig.H"

#include "gmshViewBase.H"

#if !WITH_FVSPATCHFIELD
#define fvsPatchField fvPatchField
#endif

namespace Foam
{

class gmshViews;

template <class T1, class T2>
class gmshViewSurfaceMesh
: public gmshViewBase
{

    const static label geoOffset_;

    // Reference to the mesh
    const fvMesh& mesh_;

    void writeFaceValues(char *values, const DynamicList<label>& shapeFaces,
    const GeometricField<T2, fvsPatchField, surfaceMesh>& field,
    const label timeI) const;
    void writePatchValues(char *values,
    const List<DynamicList<label> >& shapeFaces,
    const GeometricField<T2, fvsPatchField, surfaceMesh>& field,
    const label timeI, const label faceOffset) const;

protected:
    // Protected data
    const static label fieldType_;
    const static enum meshTypes meshType_;

    // Protected member functions
    void writeFacePoints(char *p, const DynamicList<label>& shapeFaces) const;
    void writePatchPoints(char *p, const List<DynamicList<label> >& shapeFaces,
    const label faceOffset) const;

public:
    // Constructor
    gmshViewSurfaceMesh(const gmshViews& views, const word& fieldName,
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
