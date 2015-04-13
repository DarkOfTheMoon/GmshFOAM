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
    gmshViewBase

Description
    The abstract base class for reading an OpenFOAM field and converts
    to Gmsh Views.

SourceFiles
    gmshViewBase.C

\*---------------------------------------------------------------------------*/

#ifndef gmshViewBase_H
#define gmshViewBase_H

#include "Time.H"
#include "fvMesh.H"

namespace Foam
{

class gmshViews;

class gmshViewBase
{
public:
    enum geoElem
    {
        nGeoTypes = 45,
        SP = 0, VP = 1, TP = 2, SL = 3, VL = 4, TL = 5,
        ST = 6, VT = 7, TT = 8, SQ = 9, VQ = 10, TQ = 11,
        SS = 12, VS = 13, TS = 14, SH = 15, VH = 16, TH = 17,
        SI = 18, VI = 19, TI = 20, SY = 21, VY = 22, TY = 23,
        SL2 = 24, VL2 = 25, TL2 = 26, ST2 = 27, VT2 = 28, TT2 = 29,
        SQ2 = 30, VQ2 = 31, TQ2 = 32, SS2 = 33, VS2 = 34, TS2 = 35,
        SH2 = 36, VH2 = 37, TH2 = 38, SI2 = 39, VI2 = 40, TI2 = 41,
        SY2 = 42, VY2 = 43, TY2 = 44
    };

    enum meshTypes { typeVolMesh, typeSurfaceMesh, typePointMesh };

    // $PostFormat
    struct gmshPostFormat
    {
        double version;
        int format, size;
    };

    // $View header
    struct gmshViewHeader
    {
        char name[256];
        int NbTimeStep;
        int Nb[nGeoTypes];
        int NbT2, t2l, NbT3, t3l;
    };

    // $View data
    struct gmshViewBuffer
    {
        char *buf[nGeoTypes];
        char *T2D, *T2C, *T3D, *T3C;
    };

protected:
    // Static data
    static label nVertices_[];
    static label nComponents_[];

    // Reference to the views
    const gmshViews& views_;

    const word fieldName_;
    Time& runTime_;

    const label verbosity_;

    // Private member functions
    inline const scalar& component(const scalar& s, const label componentI)
        const { return s; }
    inline const scalar& component(const vector& v, const label componentI)
        const { return v[componentI]; }
    inline const scalar& component(const tensor& t, const label componentI)
        const { return t[componentI]; }

    inline const scalar& component(const symmTensor& t, const label componentI)
        const
    {
        static label tab[] = {0, 1, 2, 1, 3, 4, 2, 4, 5};
        return t[tab[componentI]];
    }

public:
    // Constructor
    gmshViewBase(const gmshViews& views, const word& fieldName, Time& runTime,
    const label verbosity)
        : views_(views), fieldName_(fieldName), runTime_(runTime),
          verbosity_(verbosity)
        {}

    // Destructor - required!
    virtual ~gmshViewBase() {}

    // Pure virtual functions
    virtual const enum meshTypes meshType() const = 0;
    virtual void getPostFormat(gmshPostFormat& pF) const = 0;
    virtual void getTimeStepValues(char *tB) const = 0;
    virtual void getViewHeader(gmshViewHeader& vH) const = 0;
    virtual void getViewData(gmshViewBuffer& vB) const = 0;

    // Member functions
    void getViewSize(const gmshViewHeader& vH, label vSize[]) const;
    void write(const fileName& name, const gmshPostFormat& pF,
    const gmshViewHeader& vH, const char *tB, const gmshViewBuffer &vB) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
