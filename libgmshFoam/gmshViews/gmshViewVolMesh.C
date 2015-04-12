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
    Reads an OpenFOAM volMesh field and converts to a Gmsh View.

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "fvPatchField.H"

#include "gmshMessageStream.H"
#include "gmshViewVolMesh.H"
#include "gmshViews.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Static data * * * * * * * * * * * * * * * //

template <class T1, class T2>
const label gmshViewVolMesh<T1, T2>
::fieldType_ = (sizeof(T2) / sizeof(scalar) == 6
    ? 9 : sizeof(T2) / sizeof(scalar)); // treat a SymmTensor as a Tensor

template <class T1, class T2>
const label gmshViewVolMesh<T1, T2>
::geoOffset_ = (gmshViewVolMesh<T1, T2>::fieldType_ == 1 ? 0
: (gmshViewVolMesh<T1, T2>::fieldType_ == 3 ? 1 : 2));

template <class T1, class T2>
const gmshViewBase::meshTypes gmshViewVolMesh<T1, T2>
::meshType_ = gmshViewBase::typeVolMesh;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class T1, class T2>
gmshViewVolMesh<T1, T2>
::gmshViewVolMesh(const gmshViews& views, const word& fieldName, Time& runTime,
const fvMesh& mesh, const label verbosity)
    : gmshViewBase(views, fieldName, runTime, verbosity), mesh_(mesh)
{
    gInfo(verbosity_ >= 3) << endl << "Reading " << fieldName << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T1, class T2>
void gmshViewVolMesh<T1, T2>
::getPostFormat(gmshViewBase::gmshPostFormat& pF) const
{
    pF.version = 1.4;
    pF.format = 1; // 0: ASCII, 1: binary
    pF.size = sizeof(T1);
}

template <class T1, class T2>
void gmshViewVolMesh<T1, T2>
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
void gmshViewVolMesh<T1, T2>
::getViewHeader(gmshViewBase::gmshViewHeader& vH) const
{
    strncpy(vH.name, fieldName_.c_str(), 255);
    vH.name[255] = '\0';
    vH.NbTimeStep = views_.nTimeSteps();

    for(label typeI = 0; typeI < nGeoTypes; typeI++)
    {
        vH.Nb[typeI] = 0;
    }

    label nPatchTris = 0;
    forAll(views_.patchTris(), patchI)
    {
        nPatchTris += views_.patchTris()[patchI].size();
    }
    label nPatchQuads = 0;
    forAll(views_.patchQuads(), patchI)
    {
        nPatchQuads += views_.patchQuads()[patchI].size();
    }

    vH.Nb[ST + geoOffset_] = nPatchTris;
    vH.Nb[SQ + geoOffset_] = nPatchQuads;

    vH.Nb[SS + geoOffset_] = views_.tets().size();
    vH.Nb[SH + geoOffset_] = views_.hexes().size();
    vH.Nb[SI + geoOffset_] = views_.prisms().size();
    vH.Nb[SY + geoOffset_] = views_.pyrs().size();

    vH.NbT2 = 0; vH.t2l = 0; vH.NbT3 = 0; vH.t3l = 0;
}

template <class T1, class T2>
void gmshViewVolMesh<T1, T2>
::writeCellPoints(char *p, const labelList& shapeCells) const
{
    forAll(shapeCells, cellI)
    {
        const cellShape& shapeCell = mesh_.cellShapes()[shapeCells[cellI]];
        const pointField cellPoints(shapeCell.points(mesh_.points()));
        label pI = cellI * cellPoints.size()
            * (views_.nTimeSteps() * fieldType_ + 3);

        forAll(cellPoints, pointI)
        {
            reinterpret_cast<T1*>(p)[pI] = cellPoints[pointI].x();
            pI++;
        }
        forAll(cellPoints, pointI)
        {
            reinterpret_cast<T1*>(p)[pI] = cellPoints[pointI].y();
            pI++;
        }
        forAll(cellPoints, pointI)
        {
            reinterpret_cast<T1*>(p)[pI] = cellPoints[pointI].z();
            pI++;
        }
    }
}

template <class T1, class T2>
void gmshViewVolMesh<T1, T2>
::writePatchPoints(char *p, const List<DynamicList<label> >& shapeFaces,
const label faceOffset) const
{
    label faceI = faceOffset;
    forAll(shapeFaces, patchI)
    {
        forAll(shapeFaces[patchI], faceJ)
        {
            const face& f = mesh_.faces()[mesh_.boundaryMesh()[patchI].start()
            + shapeFaces[patchI][faceJ]];
            label pI = faceI * f.size()
                * (views_.nTimeSteps() * fieldType_ + 3);
        
            forAll(f, pointI)
            {
                reinterpret_cast<T1*>(p)[pI] = mesh_.points()[f[pointI]].x();
                pI++;
            }
            forAll(f, pointI)
            {
                reinterpret_cast<T1*>(p)[pI] = mesh_.points()[f[pointI]].y();
                pI++;
            }
            forAll(f, pointI)
            {
                reinterpret_cast<T1*>(p)[pI] = mesh_.points()[f[pointI]].z();
                pI++;
            }
            faceI++;
        }
    }
}

template <class T1, class T2>
void gmshViewVolMesh<T1, T2>
::writeCellValues(char *values, const labelList& shapeCells,
const GeometricField<T2, fvPatchField, volMesh>& field, const label timeI)
    const
{
    forAll(shapeCells, cellI)
    {
        label nP = mesh_.cellShapes()[shapeCells[cellI]].nPoints();
        label idx1 = cellI * nP * (views_.nTimeSteps() * fieldType_ + 3)
            + nP * (timeI * fieldType_ + 3);

        for(label pointI = 0; pointI < nP; pointI++)
        {
            const label idx2 = idx1 + pointI * fieldType_;

            for(label componentI = 0; componentI < fieldType_; componentI++)
            {
                reinterpret_cast<T1*>(values)[idx2 + componentI]
                    =  component(field[shapeCells[cellI]], componentI);
            }
        }
    }
}

template <class T1, class T2>
void gmshViewVolMesh<T1, T2>
::writePatchValues(char *values, const List<DynamicList<label> >& shapeFaces,
const GeometricField<T2, fvPatchField, volMesh>& field, const label timeI,
const label faceOffset)
    const
{
    label faceI = faceOffset;

    forAll(shapeFaces, patchI)
    {
        forAll(shapeFaces[patchI], faceJ)
        {
            label faceK = shapeFaces[patchI][faceJ];
            label nP = mesh_.boundaryMesh()[patchI][faceK].size();
            label idx1 = faceI * nP * (views_.nTimeSteps() * fieldType_ + 3)
                + nP * (timeI * fieldType_ + 3);

            for(label pointI = 0; pointI < nP; pointI++)
            {
                const label idx2 = idx1 + pointI * fieldType_;

                for(label componentI = 0; componentI < fieldType_;
                    componentI++)
                {
                    reinterpret_cast<T1*>(values)[idx2 + componentI] =
                        component(field.boundaryField()[patchI][faceK],
                        componentI);
                }
            }
            faceI++;
        }
    }
}

template <class T1, class T2>
void gmshViewVolMesh<T1, T2>
::getViewData(gmshViewBase::gmshViewBuffer& vB) const
{
    writeCellPoints(vB.buf[SS + geoOffset_], views_.tets());
    writeCellPoints(vB.buf[SH + geoOffset_], views_.hexes());
    writeCellPoints(vB.buf[SI + geoOffset_], views_.prisms());
    writeCellPoints(vB.buf[SY + geoOffset_], views_.pyrs());

    writePatchPoints(vB.buf[ST + geoOffset_], views_.patchTris(), 0);
    writePatchPoints(vB.buf[SQ + geoOffset_], views_.patchQuads(), 0);

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
        GeometricField<T2, fvPatchField, volMesh>* fieldPtr = NULL;
        try
        {
            fieldObjectPtr = new IOobject(fieldName_, runTime_.timeName(),
            mesh_, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE);
            fieldPtr = new GeometricField<T2, fvPatchField, volMesh>
                (*fieldObjectPtr, mesh_,
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
        GeometricField<T2, fvPatchField, volMesh>& field = *fieldPtr;

        if(!fieldObject.headerOk())
        {
            field = pTraits<T2>::zero;
            gWarning(verbosity_ >= 1) << "Field " << fieldName_
                << " not found at t = " << runTime_.timeName()
                << ". Treating the field values as all zero." << endl;
        }

        const label timeI = curTime - views_.startTime();

        writeCellValues(vB.buf[SS + geoOffset_], views_.tets(), field, timeI);
        writeCellValues(vB.buf[SH + geoOffset_], views_.hexes(), field, timeI);
        writeCellValues(vB.buf[SI + geoOffset_], views_.prisms(), field,
        timeI);
        writeCellValues(vB.buf[SY + geoOffset_], views_.pyrs(), field, timeI);

        writePatchValues(vB.buf[ST + geoOffset_], views_.patchTris(), field,
        timeI, 0);
        writePatchValues(vB.buf[SQ + geoOffset_], views_.patchQuads(), field,
        timeI, 0);

        delete fieldPtr;
        delete fieldObjectPtr;
    }
}

// * * * * * * * * * * * * * * * * Instantiators * * * * * * * * * * * * * * //

// double precision Gmsh view / volScalarField
template class gmshViewVolMesh<double, scalar>;
// double precision Gmsh view / volVectorField
template class gmshViewVolMesh<double, vector>;
#if WITH_SYMMTENSOR
// double precision Gmsh view / volSymmTensorField
template class gmshViewVolMesh<double, symmTensor>;
#endif
// double precision Gmsh view / volTensorField
template class gmshViewVolMesh<double, tensor>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
