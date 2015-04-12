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
    Reads an OpenFOAM volMesh/surfaceMesh field and converts to a Gmsh View.

\*---------------------------------------------------------------------------*/

#include "OFstream.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "fvPatchField.H"

#include "gmshMessageStream.H"
#include "gmshView.H"
#include "gmshViews.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Static data * * * * * * * * * * * * * * * //

template <class T1, class T2, class T4>
const label gmshView<T1, T2, T4>
::fieldType_ = sizeof(T2) / sizeof(scalar);

template <class T1, class T2, class T4>
const label gmshView<T1, T2, T4>
::geoOffset_ = (gmshView<T1, T2, T4>::fieldType_ == 1 ? 0
: (gmshView<T1, T2, T4>::fieldType_ == 3 ? 1 : 2));

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class T1, class T2, class T4>
gmshView<T1, T2, T4>
::gmshView(const gmshViews& views, const word& fieldName, Time& runTime,
const fvMesh& mesh, const label verbosity)
    : gmshViewBase(views, fieldName, runTime, verbosity), mesh_(mesh)
{

    const word& typeName = GeometricField<T2, fvPatchField, T4>::typeName;

    if(typeName == "volScalarField" || typeName == "volVectorField"
    || typeName == "volTensorField")
    {
        meshType_ = typeVolMesh;
    }
    else if(typeName == "surfaceScalarField"
    || typeName == "surfaceVectorField" || typeName == "surfaceTensorField")
    {
        meshType_ = typeSurfaceMesh;
    }
    else
    {
        // usually this should not happen
        gSeriousError(verbosity_ >= 1) << "Unhandled field type; type = "
            << typeName << endl;
    }

    gInfo(verbosity_ >= 3) << endl << "Reading " << fieldName << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T1, class T2, class T4>
void gmshView<T1, T2, T4>
::getPostFormat(gmshViewBase::gmshPostFormat& pF) const
{
    pF.version = 1.4;
    pF.format = 1; // 0: ASCII, 1: binary
    pF.size = sizeof(T1);
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

    if(meshType_ == typeVolMesh)
    {
        vH.Nb[ST + geoOffset_] = nPatchTris;
        vH.Nb[SQ + geoOffset_] = nPatchQuads;

        vH.Nb[SS + geoOffset_] = views_.tets().size();
        vH.Nb[SH + geoOffset_] = views_.hexes().size();
        vH.Nb[SI + geoOffset_] = views_.prisms().size();
        vH.Nb[SY + geoOffset_] = views_.pyrs().size();
    }
    else if(meshType_ == typeSurfaceMesh)
    {
        vH.Nb[ST + geoOffset_] = views_.tris().size() + nPatchTris;
        vH.Nb[SQ + geoOffset_] = views_.quads().size() + nPatchQuads;
    }

    vH.NbT2 = 0; vH.t2l = 0; vH.NbT3 = 0; vH.t3l = 0;
}

template <class T1, class T2, class T4>
void gmshView<T1, T2, T4>
::getTimeStepValues(char *tB) const
{
    for(label curTime = views_.startTime(); curTime <= views_.endTime();
        curTime++)
    {
        reinterpret_cast<T1*>(tB)[curTime - views_.startTime()]
            = views_.timeList()[curTime].value();
    }
}

template <class T1, class T2, class T4>
void gmshView<T1, T2, T4>
::writeCellPoints(char *p, const labelList& shapeCells) const
{
    if(meshType_ != typeVolMesh)
    {
        return;
    }

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

template <class T1, class T2, class T4>
void gmshView<T1, T2, T4>
::writeFacePoints(char *p, const DynamicList<label>& shapeFaces) const
{
    if(meshType_ != typeSurfaceMesh)
    {
        return;
    }

    forAll(shapeFaces, faceI)
    {
        const face& f = mesh_.faces()[shapeFaces[faceI]];
        label pI = faceI * f.size() * (views_.nTimeSteps() * fieldType_ + 3);
        
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
    }
}

template <class T1, class T2, class T4>
void gmshView<T1, T2, T4>
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

template <class T1, class T2, class T4>
void gmshView<T1, T2, T4>
::writeCellValues(char *values, const labelList& shapeCells,
const GeometricField<T2, fvPatchField, T4>& field, const label timeI)
    const
{
    if(meshType_ != typeVolMesh)
    {
        return;
    }

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

template <class T1, class T2, class T4>
void gmshView<T1, T2, T4>
::writeFaceValues(char *values, const DynamicList<label>& shapeFaces,
const GeometricField<T2, fvPatchField, T4>& field, const label timeI)
    const
{
    if(meshType_ != typeSurfaceMesh)
    {
        return;
    }

    forAll(shapeFaces, faceI)
    {
        label nP = mesh_.faces()[shapeFaces[faceI]].size();
        label idx1 = faceI * nP * (views_.nTimeSteps() * fieldType_ + 3)
            + nP * (timeI * fieldType_ + 3);

        for(label pointI = 0; pointI < nP; pointI++)
        {
            const label idx2 = idx1 + pointI * fieldType_;

            for(label componentI = 0; componentI < fieldType_; componentI++)
            {
                reinterpret_cast<T1*>(values)[idx2 + componentI]
                    =  component(field[shapeFaces[faceI]], componentI);
            }
        }
    }
}

template <class T1, class T2, class T4>
void gmshView<T1, T2, T4>
::writePatchValues(char *values, const List<DynamicList<label> >& shapeFaces,
const GeometricField<T2, fvPatchField, T4>& field, const label timeI,
const label faceOffset)
    const
{
    if(meshType_ != typeVolMesh && meshType_ != typeSurfaceMesh)
    {
        return;
    }

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

template <>
void gmshView<double, vector, pointMesh>
::writePointValues(char *values,
const GeometricField<vector, pointPatchField, pointMesh>& field,
const label timeI)
    const
{
    if(meshType_ != typePointMesh)
    {
        return;
    }

    forAll(mesh_.points(), pointI)
    {
        const label idx = pointI * (views_.nTimeSteps() * fieldType_ + 3)
            + (timeI * fieldType_ + 3);

        for(label componentI = 0; componentI < fieldType_; componentI++)
        {
            reinterpret_cast<double*>(values)[idx + componentI]
                =  component(field[pointI], componentI);
        }
    }
    // currently no support for patch values
}

template <class T1, class T2, class T4>
void gmshView<T1, T2, T4>
::getViewData(gmshViewBase::gmshViewBuffer& vB) const
{
    writeFacePoints(vB.buf[ST + geoOffset_], views_.tris());
    writeFacePoints(vB.buf[SQ + geoOffset_], views_.quads());

    writeCellPoints(vB.buf[SS + geoOffset_], views_.tets());
    writeCellPoints(vB.buf[SH + geoOffset_], views_.hexes());
    writeCellPoints(vB.buf[SI + geoOffset_], views_.prisms());
    writeCellPoints(vB.buf[SY + geoOffset_], views_.pyrs());

    writePatchPoints(vB.buf[ST + geoOffset_], views_.patchTris(),
    (meshType_ == typeSurfaceMesh ? views_.tris().size() : 0));
    writePatchPoints(vB.buf[SQ + geoOffset_], views_.patchQuads(),
    (meshType_ == typeSurfaceMesh ? views_.quads().size() : 0));

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
        IOobject::READ_IF_PRESENT, IOobject::NO_WRITE);
        GeometricField<T2, fvPatchField, T4> field(fieldObject, mesh_,
        dimensionSet(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

        if(!fieldObject.headerOk())
        {
            field = pTraits<T2>::zero;
            gSeriousError(verbosity_ >= 1) << "Field " << fieldName_
                << " not found at t = " << runTime_.timeName()
                << ". Treating the field values as all zero." << endl;
        }

        const label timeI = curTime - views_.startTime();

        writeFaceValues(vB.buf[ST + geoOffset_], views_.tris(), field, timeI);
        writeFaceValues(vB.buf[SQ + geoOffset_], views_.quads(), field, timeI);

        writeCellValues(vB.buf[SS + geoOffset_], views_.tets(), field, timeI);
        writeCellValues(vB.buf[SH + geoOffset_], views_.hexes(), field, timeI);
        writeCellValues(vB.buf[SI + geoOffset_], views_.prisms(), field,
        timeI);
        writeCellValues(vB.buf[SY + geoOffset_], views_.pyrs(), field, timeI);

        writePatchValues(vB.buf[ST + geoOffset_], views_.patchTris(), field,
        timeI, (meshType_ == gmshViewBase::typeSurfaceMesh
        ? views_.tris().size() : 0));
        writePatchValues(vB.buf[SQ + geoOffset_], views_.patchQuads(), field,
        timeI, (meshType_ == gmshViewBase::typeSurfaceMesh
        ? views_.quads().size() : 0));
    }
}

template <class T1, class T2, class T4>
void gmshView<T1, T2, T4>
::write(const fileName& name, const gmshPostFormat& pF,
const gmshViewHeader& vH, const char* tB, const gmshViewBuffer &vB) const
{
    OFstream vFile(name, OFstream::BINARY);
    OStringStream oStr;

    gInfo(verbosity_ >= 3) << "    Writing view" << endl;
    gInfo(verbosity_ >= 4) << "        Writing header" << endl;

    // write header
    oStr << "$PostFormat" << endl
        << "1.4 0 " << pF.size << endl
        << "$EndPostFormat" << endl;

    oStr << "$View" << endl;
    oStr << vH.name << " " << vH.NbTimeStep << endl;
    for(label typeI = 0; typeI < gmshViewBase::nGeoTypes / 3; typeI++)
    {
        for(label typeJ = 0; typeJ < 3; typeJ++)
        {
            oStr << vH.Nb[typeI * 3 + typeJ] << " ";
        }
        oStr << endl;
    }
    oStr << vH.NbT2 << " " << vH.t2l << " "
        << vH.NbT3 << " " << vH.t3l << endl;

    vFile << oStr.str().c_str();
    for(label timeI = 0; timeI < vH.NbTimeStep; timeI++)
    {
        vFile << reinterpret_cast<const T1*>(tB)[timeI] << " ";
    }

    gInfo(verbosity_ >= 4) << "        Writing elements" << endl;

    vFile << endl;
    for(label typeI = 0; typeI < gmshViewBase::nGeoTypes; typeI++)
    {
        label elemI = 0;
        for(label cellI = 0; cellI < vH.Nb[typeI]; cellI++)
        {
            for(label componentI = 0; componentI < 3; componentI++)
            {
                for(label pointI = 0; pointI < nVertices_[typeI / 3]; pointI++)
                {
                    vFile << reinterpret_cast<T1*>(vB.buf[typeI])[elemI]
                        << " ";
                    elemI++;
                }
                vFile << endl;
            }
            for(label timeI = 0; timeI < vH.NbTimeStep; timeI++)
            {
                for(label pointI = 0; pointI < nVertices_[typeI / 3]; pointI++)
                {
                    for(label componentI = 0;
                        componentI < nComponents_[typeI % 3]; componentI++)
                    {
                        vFile << reinterpret_cast<T1*>(vB.buf[typeI])[elemI]
                            << " ";
                        elemI++;
                    }
                    vFile << endl;
                }
            }
        }
    }

    // 2D texts
    for(label textI = 0, elemI = 0; textI < vH.NbT2; textI++)
    {
        for(label numI = 0; numI < 4; numI++)
        {
            vFile << reinterpret_cast<T1*>(vB.T2D)[elemI] << " ";
            elemI++;
        }
        vFile << endl;
    }
    for(label charI = 0; charI < vH.t2l; charI++)
    {
        vFile << vB.T2C[charI];
    }
    vFile << endl;

    // 3D texts
    for(label textI = 0, elemI = 0; textI < vH.NbT3; textI++)
    {
        for(label numI = 0; numI < 5; numI++)
        {
            vFile << reinterpret_cast<T1*>(vB.T3D)[elemI] << " ";
            elemI++;
        }
        vFile << endl;
    }
    for(label charI = 0; charI < vH.t3l; charI++)
    {
        vFile << vB.T3C[charI];
    }
    vFile << endl;

    vFile << "$EndView" << endl;

    gInfo(verbosity_ >= 4) << "        ... done." << endl;
}

// * * * * * * * * * * * * * * * * Instantiators * * * * * * * * * * * * * * //

// double precision Gmsh view / volScalarField
template class gmshView<double, scalar, volMesh>;
// double precision Gmsh view / volVectorField
template class gmshView<double, vector, volMesh>;
// double precision Gmsh view / volVectorField
template class gmshView<double, tensor, volMesh>;

// double precision Gmsh view / surfaceScalarField
template class gmshView<double, scalar, surfaceMesh>;
// double precision Gmsh view / surfaceVectorField
template class gmshView<double, vector, surfaceMesh>;
// double precision Gmsh view / surfaceVectorField
template class gmshView<double, tensor, surfaceMesh>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
