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
    Reads an OpenFOAM field and converts to a Gmsh View.

\*---------------------------------------------------------------------------*/

#include "OFstream.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "fvPatchField.H"

#include "gmshMessageStream.H"
#include "gmshViewMeshMotion.H"
#include "gmshViews.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class T1, class T4>
gmshViewMeshMotion<T1, T4>
::gmshViewMeshMotion(const gmshViews& views, const word& fieldName,
Time& runTime, const fvMesh& mesh, const label verbosity)
    :gmshView<T1, vector, T4>(views, fieldName, runTime, mesh, verbosity),
     mesh_(mesh), views_(views),
     fieldType_(gmshView<T1, vector, T4>::fieldType_),
     meshType_(gmshView<T1, vector, T4>::meshType_), runTime_(runTime),
     verbosity_(verbosity)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T1, class T4>
void gmshViewMeshMotion<T1, T4>
::writeFaceMotions(char *values, const DynamicList<label>& shapeFaces,
pointField& meshMotion, const label timeI)
    const
{
    if(meshType_ != gmshViewBase::typeSurfaceMesh || fieldType_ != 3)
    {
        return;
    }

    forAll(shapeFaces, faceI)
    {
        const face& f = mesh_.faces()[shapeFaces[faceI]];
        label nP = f.size();
        label idx1 = faceI * nP * (views_.nTimeSteps() * fieldType_ + 3)
            + nP * (timeI * fieldType_ + 3);

        for(label pointI = 0; pointI < nP; pointI++)
        {
            const label idx2 = idx1 + pointI * fieldType_;

            for(label componentI = 0; componentI < fieldType_; componentI++)
            {
                reinterpret_cast<T1*>(values)[idx2 + componentI]
                    =  component(meshMotion[f[pointI]], componentI);
            }
        }
    }
}

template <class T1, class T4>
void gmshViewMeshMotion<T1, T4>
::writeCellMotions(char *values, const labelList& shapeCells,
pointField& meshMotion, const label timeI) const
{
    if(meshType_ != gmshViewBase::typeVolMesh || fieldType_ != 3)
    {
        return;
    }

    forAll(shapeCells, cellI)
    {
        const cellShape& c = mesh_.cellShapes()[shapeCells[cellI]];
        label nP = c.size();
        label idx1 = cellI * nP * (views_.nTimeSteps() * fieldType_ + 3)
            + nP * (timeI * fieldType_ + 3);

        for(label pointI = 0; pointI < nP; pointI++)
        {
            const label idx2 = idx1 + pointI * fieldType_;

            for(label componentI = 0; componentI < fieldType_; componentI++)
            {
                reinterpret_cast<T1*>(values)[idx2 + componentI]
                    =  component(meshMotion[c[pointI]], componentI);
            }
        }
    }
}

template <class T1, class T4>
void gmshViewMeshMotion<T1, T4>
::writePatchMotions(char *values, const List<DynamicList<label> >& shapeFaces,
const pointField& meshMotion, const label timeI, const label faceOffset)
    const
{
    if((meshType_ != gmshViewBase::typeVolMesh
       && meshType_ != gmshViewBase::typeSurfaceMesh)
    || fieldType_ != 3)
    {
        return;
    }

    label faceI = faceOffset;

    forAll(shapeFaces, patchI)
    {
        forAll(shapeFaces[patchI], faceJ)
        {
            label faceK = shapeFaces[patchI][faceJ];
            const face& f = mesh_.boundaryMesh()[patchI][faceK];
            label nP = f.size();
            label idx1 = faceI * nP * (views_.nTimeSteps() * fieldType_ + 3)
                + nP * (timeI * fieldType_ + 3);

            for(label pointI = 0; pointI < nP; pointI++)
            {
                const label idx2 = idx1 + pointI * fieldType_;

                for(label componentI = 0; componentI < fieldType_;
                    componentI++)
                {
                    reinterpret_cast<T1*>(values)[idx2 + componentI] =
                        component(meshMotion[f[pointI]], componentI);
                }
            }
            faceI++;
        }
    }
}

template <class T1, class T4>
void gmshViewMeshMotion<T1, T4>
::getViewData(gmshViewBase::gmshViewBuffer& vB) const
{
    if(fieldType_ != 3)
    {
        gSeriousError(verbosity_ >= 1) << "unsupported template instance for"
            "handling mesh motions. fieldType = " << fieldType_ << endl;
        return;
    }

    pointField meshMotion(mesh_.points().size(), vector::zero);

    gmshView<T1, vector, T4>::
        writeFacePoints(vB.buf[gmshViewBase::VT], views_.tris());
    gmshView<T1, vector, T4>::
        writeFacePoints(vB.buf[gmshViewBase::VQ], views_.quads());

    gmshView<T1, vector, T4>::
        writeCellPoints(vB.buf[gmshViewBase::VS], views_.tets());
    gmshView<T1, vector, T4>::
        writeCellPoints(vB.buf[gmshViewBase::VH], views_.hexes());
    gmshView<T1, vector, T4>::
        writeCellPoints(vB.buf[gmshViewBase::VI], views_.prisms());
    gmshView<T1, vector, T4>::
        writeCellPoints(vB.buf[gmshViewBase::VY], views_.pyrs());

    gmshView<T1, vector, T4>::
        writePatchPoints(vB.buf[gmshViewBase::VT], views_.patchTris(),
        (meshType_ == gmshViewBase::typeSurfaceMesh
        ? views_.tris().size() : 0));
    gmshView<T1, vector, T4>::
        writePatchPoints(vB.buf[gmshViewBase::VQ], views_.patchQuads(),
        (meshType_ == gmshViewBase::typeSurfaceMesh
        ? views_.quads().size() : 0));

    for(label curTime = views_.startTime(); curTime <= views_.endTime();
        curTime++)
    {
        gInfo(verbosity_ >= 4) << "    Reading time = "
            << views_.timeList()[curTime].name();
        gInfo(verbosity_ >= 5) << " (" << curTime - views_.startTime() + 1
            << "/" << views_.nTimeSteps() << ")";
        gInfo(verbosity_ >= 4) << endl;

        runTime_.setTime(views_.timeList()[curTime], curTime);
        IOobject ioPoints("points", runTime_.timeName(), polyMesh::typeName,
        mesh_, IOobject::MUST_READ, IOobject::NO_WRITE);
        if(ioPoints.headerOk())
        {
            gInfo(verbosity_ >= 4)
                << "        Applying mesh movement at t = "
                << views_.timeList()[curTime].name() << endl;

            pointIOField newPoints(ioPoints);
            meshMotion = newPoints - mesh_.points();
        }

        const label timeI = curTime - views_.startTime();

        writeFaceMotions(vB.buf[gmshViewBase::VT], views_.tris(), meshMotion,
        timeI);
        writeFaceMotions(vB.buf[gmshViewBase::VQ], views_.quads(), meshMotion,
        timeI);

        writeCellMotions(vB.buf[gmshViewBase::VS], views_.tets(), meshMotion,
        timeI);
        writeCellMotions(vB.buf[gmshViewBase::VH], views_.hexes(), meshMotion,
        timeI);
        writeCellMotions(vB.buf[gmshViewBase::VI], views_.prisms(), meshMotion,
        timeI);
        writeCellMotions(vB.buf[gmshViewBase::VY], views_.pyrs(), meshMotion,
        timeI);

        writePatchMotions(vB.buf[gmshViewBase::VT], views_.patchTris(),
        meshMotion, timeI, (meshType_ == gmshViewBase::typeSurfaceMesh
        ? views_.tris().size() : 0));
        writePatchMotions(vB.buf[gmshViewBase::VQ], views_.patchQuads(),
        meshMotion, timeI, (meshType_ == gmshViewBase::typeSurfaceMesh
        ? views_.quads().size() : 0));
    }
}

// * * * * * * * * * * * * * * * * Instantiators * * * * * * * * * * * * * * //

// double precision Gmsh view / volMesh
template class gmshViewMeshMotion<double, volMesh>;
// double precision Gmsh view / surfaceMesh
template class gmshViewMeshMotion<double, surfaceMesh>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
