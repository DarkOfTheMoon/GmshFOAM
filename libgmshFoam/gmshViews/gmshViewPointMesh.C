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
    Reads an OpenFOAM pointMesh field and converts to a Gmsh View.

    The coding is still unfinished. This file is considered as a scratch file.

\*---------------------------------------------------------------------------*/

#include "pointFields.H"

template <>
gmshView<double, vector, pointMesh>
::gmshView(const gmshViews& views, const word& fieldName, Time& runTime,
const fvMesh& mesh, const label verbosity)
    : gmshViewBase(views, fieldName, runTime, verbosity), mesh_(mesh)
{

    const word& typeName
        = GeometricField<double, pointPatchField, pointMesh>::typeName;

    if(typeName == "pointScalarField" || typeName == "pointVectorField"
    || typeName == "pointTensorField")
    {
        meshType_ = typePointMesh;
    }
    else
    {
        // usually this should not happen
        gSeriousError(verbosity_ >= 1) << "Unhandled field type; type = "
            << typeName << endl;
    }

    gInfo(verbosity_ >= 3) << endl << "Reading " << fieldName << endl;
}

template <class T1, class T2, class T4>
void gmshView<T1, T2, T4>
::getViewHeader(gmshViewBase::gmshViewHeader& vH) const
{
    strncpy(vH.name, fieldName_.c_str(), 255);
    vH.name[255] = '\0';
    vH.NbTimeStep = views_.nTimeSteps();

    for(label typeI = 0; typeI < nGeoTypes; typeI++)
    {
        vH.Nb[typeI] = 0;
    }

    const static label type = sizeof(T2) / sizeof(scalar);
    const static label offset = (type == 1 ? 0 : (type == 3 ? 1 : 2));

    if(meshType_ == typePointMesh)
    {
        vH.Nb[SP + offset] = mesh_.points().size();
        //vH.Nb[ST + offset] = nPatchTris;
        //vH.Nb[SQ + offset] = nPatchQuads;
    }

    vH.NbT2 = 0; vH.t2l = 0; vH.NbT3 = 0; vH.t3l = 0;
}

template <class T1, class T2, class T4>
void gmshView<T1, T2, T4>
::writePointPoints(char *p) const
{
    if(meshType_ != typePointMesh)
    {
        return;
    }

    forAll(mesh_.points(), pointI)
    {
        const label pI = pointI * (views_.nTimeSteps() * fieldType_ + 3);

        reinterpret_cast<T1*>(p)[pI] = mesh_.points()[pointI].x();
        reinterpret_cast<T1*>(p)[pI + 1] = mesh_.points()[pointI].y();
        reinterpret_cast<T1*>(p)[pI + 2] = mesh_.points()[pointI].z();
    }
}

template <>
void gmshView<double, vector, pointMesh>
::getViewData(gmshViewBase::gmshViewBuffer& vB) const
{
    writePointPoints(vB.buf[VP]);

    for(label curTime = views_.startTime(); curTime <= views_.endTime();
        curTime++)
    {
        gInfo(verbosity_ >= 4) << "    Reading time = "
            << views_.timeList()[curTime].name();
        gInfo(verbosity_ >= 5) << " (" << curTime - views_.startTime() + 1
            << "/" << views_.nTimeSteps() << ")";
        gInfo(verbosity_ >= 4) << endl;

        runTime_.setTime(views_.timeList()[curTime], curTime);

        IOobject fieldObject(fieldName_, runTime_.timeName(), mesh_,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE);

        pointMesh pMesh(mesh_);
        GeometricField<vector, pointPatchField, pointMesh>
            field(fieldObject, pMesh,
            dimensionSet(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        if(!fieldObject.headerOk())
        {
            field = pTraits<T2>::zero;
            gSeriousError(verbosity_ >= 1) << "Field " << fieldName_
                << " not found at t = " << runTime_.timeName()
                << ". Treating the field values as all zero." << endl;
        }

        const label timeI = curTime - views_.startTime();

        writePointValues(vB.buf[VP], field, timeI);
    }
}

// * * * * * * * * * * * * * * * * Instantiators * * * * * * * * * * * * * * //

// double precision Gmsh view / pointVectorField
template class gmshView<double, vector, pointMesh>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
