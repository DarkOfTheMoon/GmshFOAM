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
    Reads an OpenFOAM volMesh motion and converts to a Gmsh View.

\*---------------------------------------------------------------------------*/

#include "OFstream.H"
#include "volFields.H"
#include "fvPatchField.H"

#include "gmshMessageStream.H"
#include "gmshViewVolMeshMotion.H"
#include "gmshViews.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class T1>
gmshViewVolMeshMotion<T1>
::gmshViewVolMeshMotion(const gmshViews& views, const word& fieldName,
Time& runTime, const fvMesh& mesh, const label verbosity)
    :gmshViewVolMesh<T1, vector>(views, fieldName, runTime, mesh, verbosity),
     mesh_(mesh), views_(views),
     fieldType_(gmshViewVolMesh<T1, vector>::fieldType_), runTime_(runTime),
     verbosity_(verbosity)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T1>
void gmshViewVolMeshMotion<T1>
::writeCellMotions(char *values, const labelList& shapeCells,
pointField& meshMotion, const label timeI) const
{
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

template <class T1>
void gmshViewVolMeshMotion<T1>
::writePatchMotions(char *values, const List<DynamicList<label> >& shapeFaces,
const pointField& meshMotion, const label timeI)
    const
{
    label faceI = 0;
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

template <class T1>
void gmshViewVolMeshMotion<T1>
::getViewData(gmshViewBase::gmshViewBuffer& vB) const
{
    pointField meshMotion(mesh_.points().size(), vector::zero);

    gmshViewVolMesh<T1, vector>::
        writeCellPoints(vB.buf[gmshViewBase::VS], views_.tets());
    gmshViewVolMesh<T1, vector>::
        writeCellPoints(vB.buf[gmshViewBase::VH], views_.hexes());
    gmshViewVolMesh<T1, vector>::
        writeCellPoints(vB.buf[gmshViewBase::VI], views_.prisms());
    gmshViewVolMesh<T1, vector>::
        writeCellPoints(vB.buf[gmshViewBase::VY], views_.pyrs());

    gmshViewVolMesh<T1, vector>::
        writePatchPoints(vB.buf[gmshViewBase::VT], views_.patchTris(), 0);
    gmshViewVolMesh<T1, vector>::
        writePatchPoints(vB.buf[gmshViewBase::VQ], views_.patchQuads(), 0);

    for(label curTime = views_.startTime(); curTime <= views_.endTime();
        curTime++)
    {
        gInfo(verbosity_ >= 4) << "    Reading time = "
            << views_.timeList()[curTime].name();
        gInfo(verbosity_ >= 5) << " (" << curTime - views_.startTime() + 1
            << "/" << views_.nTimeSteps() << ")";
        gInfo(verbosity_ >= 4) << endl;

        runTime_.setTime(views_.timeList()[curTime], curTime);

        try
        {
            IOobject ioPoints("points", runTime_.timeName(),
            polyMesh::typeName, mesh_, IOobject::MUST_READ,
            IOobject::NO_WRITE);

            if(ioPoints.headerOk())
            {
                gInfo(verbosity_ >= 4)
                    << "        Applying mesh movement at t = "
                    << views_.timeList()[curTime].name() << endl;

                pointIOField newPoints(ioPoints);
                meshMotion = newPoints - mesh_.points();
            }
        }
        catch(error& e)
        {
            gSeriousError(verbosity_ >= 1) << e.message().c_str() << endl;
            // treat point locations same as the ones in previous timestep
        }

        const label timeI = curTime - views_.startTime();

        writeCellMotions(vB.buf[gmshViewBase::VS], views_.tets(), meshMotion,
        timeI);
        writeCellMotions(vB.buf[gmshViewBase::VH], views_.hexes(), meshMotion,
        timeI);
        writeCellMotions(vB.buf[gmshViewBase::VI], views_.prisms(), meshMotion,
        timeI);
        writeCellMotions(vB.buf[gmshViewBase::VY], views_.pyrs(), meshMotion,
        timeI);

        writePatchMotions(vB.buf[gmshViewBase::VT], views_.patchTris(),
        meshMotion, timeI);
        writePatchMotions(vB.buf[gmshViewBase::VQ], views_.patchQuads(),
        meshMotion, timeI);
    }
}

// * * * * * * * * * * * * * * * * Instantiators * * * * * * * * * * * * * * //

// double precision Gmsh view / volMesh
template class gmshViewVolMeshMotion<double>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
