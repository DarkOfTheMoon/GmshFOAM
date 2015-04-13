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

\*---------------------------------------------------------------------------*/

#include "pointFields.H"
#include "gmshMessageStream.H"
#include "gmshViewPointMesh.H"
#include "gmshViews.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Static data * * * * * * * * * * * * * * * //

template <class T1, class T2>
const label gmshViewPointMesh<T1, T2>
::fieldType_ = (sizeof(T2) / sizeof(scalar) == 6
    ? 9 : sizeof(T2) / sizeof(scalar)); // treat a SymmTensor as a Tensor

template <class T1, class T2>
const label gmshViewPointMesh<T1, T2>
::geoOffset_ = (gmshViewPointMesh<T1, T2>::fieldType_ == 1 ? 0
: (gmshViewPointMesh<T1, T2>::fieldType_ == 3 ? 1 : 2));

template <class T1, class T2>
const gmshViewBase::meshTypes gmshViewPointMesh<T1, T2>
::meshType_ = gmshViewBase::typePointMesh;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class T1, class T2>
gmshViewPointMesh<T1, T2>
::gmshViewPointMesh(const gmshViews& views, const word& fieldName,
Time& runTime, const fvMesh& mesh, const label verbosity)
    : gmshViewBase(views, fieldName, runTime, verbosity), mesh_(mesh)
{
    gInfo(verbosity_ >= 3) << endl << "Reading " << fieldName << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T1, class T2>
void gmshViewPointMesh<T1, T2>
::getPostFormat(gmshViewBase::gmshPostFormat& pF) const
{
    pF.version = 1.4;
    pF.format = 1; // 0: ASCII, 1: binary
    pF.size = sizeof(T1);
}

template <class T1, class T2>
void gmshViewPointMesh<T1, T2>
::getTimeStepValues(char *tB) const
{
    for(label curTime = views_.startTime(); curTime <= views_.endTime();
        curTime++)
    {
        reinterpret_cast<T1*>(tB)[curTime - views_.startTime()]
            = views_.timeList()[curTime].value();
    }
}

template <class T1, class T2>
void gmshViewPointMesh<T1, T2>
::getViewHeader(gmshViewBase::gmshViewHeader& vH) const
{
    strncpy(vH.name, fieldName_.c_str(), 255);
    vH.name[255] = '\0';
    vH.NbTimeStep = views_.nTimeSteps();

    for(label typeI = 0; typeI < nGeoTypes; typeI++)
    {
        vH.Nb[typeI] = 0;
    }

    vH.Nb[SP + geoOffset_] = mesh_.points().size();

    vH.NbT2 = 0; vH.t2l = 0; vH.NbT3 = 0; vH.t3l = 0;
}

template <class T1, class T2>
void gmshViewPointMesh<T1, T2>
::writePointPoints(char *p) const
{
    forAll(mesh_.points(), pointI)
    {
        const label pI = pointI * (views_.nTimeSteps() * fieldType_ + 3);

        reinterpret_cast<T1*>(p)[pI] = mesh_.points()[pointI].x();
        reinterpret_cast<T1*>(p)[pI + 1] = mesh_.points()[pointI].y();
        reinterpret_cast<T1*>(p)[pI + 2] = mesh_.points()[pointI].z();
    }
}

template <class T1, class T2>
void gmshViewPointMesh<T1, T2>
::writePointValues(char *values,
const GeometricField<T2, pointPatchField, pointMesh>& field, const label timeI)
    const
{
    forAll(mesh_.points(), pointI)
    {
        const label idx = pointI * (views_.nTimeSteps() * fieldType_ + 3)
            + (timeI * fieldType_ + 3);
        for(label componentI = 0; componentI < fieldType_; componentI++)
        {
            reinterpret_cast<T1*>(values)[idx + componentI]
                =  component(field[pointI], componentI);
        }
    }
}

template <class T1, class T2>
void gmshViewPointMesh<T1, T2>
::getViewData(gmshViewBase::gmshViewBuffer& vB) const
{
    writePointPoints(vB.buf[SP + geoOffset_]);

    const pointMesh pMesh(mesh_);
    for(label curTime = views_.startTime(); curTime <= views_.endTime();
        curTime++)
    {
        gInfo(verbosity_ >= 4) << "    Reading time = "
            << views_.timeList()[curTime].name();
        gInfo(verbosity_ >= 5) << " (" << curTime - views_.startTime() + 1
            << "/" << views_.nTimeSteps() << ")";
        gInfo(verbosity_ >= 4) << endl;

        runTime_.setTime(views_.timeList()[curTime], curTime);

        IOobject* fieldObjectPtr = NULL;
        GeometricField<T2, pointPatchField, pointMesh>* fieldPtr = NULL;
        try
        {
            fieldObjectPtr = new IOobject(fieldName_, runTime_.timeName(),
            mesh_, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE);
            fieldPtr = new GeometricField<T2, pointPatchField, pointMesh>
            (*fieldObjectPtr, pMesh,
            dimensionSet(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        }
        catch(error& e)
        {
            delete fieldObjectPtr;
            delete fieldPtr;
            gSeriousError(verbosity_ >= 1) << e.message().c_str()
                << " at t = " << runTime_.timeName()
                << ". Treating the field values as all zero." << endl;
            continue;
        }
        IOobject& fieldObject = *fieldObjectPtr;
        GeometricField<T2, pointPatchField, pointMesh>& field = *fieldPtr;

        if(!fieldObject.headerOk())
        {
            field = pTraits<T2>::zero;
            gWarning(verbosity_ >= 1) << "Field " << fieldName_
                << " not found at t = " << runTime_.timeName()
                << ". Treating the field values as all zero." << endl;
        }

        const label timeI = curTime - views_.startTime();

        writePointValues(vB.buf[SP + geoOffset_], field, timeI);

        delete fieldPtr;
        delete fieldObjectPtr;
    }
}

// * * * * * * * * * * * * * * * * Instantiators * * * * * * * * * * * * * * //

// double precision Gmsh view / pointScalarField
template class gmshViewPointMesh<double, scalar>;
// double precision Gmsh view / pointVectorField
template class gmshViewPointMesh<double, vector>;
// double precision Gmsh view / pointSymmTensorField
template class gmshViewPointMesh<double, symmTensor>;
// double precision Gmsh view / pointTensorField
template class gmshViewPointMesh<double, tensor>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
