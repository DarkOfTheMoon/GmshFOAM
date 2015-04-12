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
    Reads an OpenFOAM pointMesh motion and converts to a Gmsh View.

\*---------------------------------------------------------------------------*/

#include "OFstream.H"
#include "pointFields.H"

#include "gmshMessageStream.H"
#include "gmshViewPointMeshMotion.H"
#include "gmshViews.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class T1>
gmshViewPointMeshMotion<T1>
::gmshViewPointMeshMotion(const gmshViews& views, const word& fieldName,
Time& runTime, const fvMesh& mesh, const label verbosity)
    :gmshViewPointMesh<T1, vector>(views, fieldName, runTime, mesh, verbosity),
     mesh_(mesh), views_(views),
     fieldType_(gmshViewPointMesh<T1, vector>::fieldType_), runTime_(runTime),
     verbosity_(verbosity)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T1>
void gmshViewPointMeshMotion<T1>
::writePointMotions(char *values, pointField& meshMotion, const label timeI)
    const
{
    forAll(mesh_.points(), pointI)
    {
        const label idx = pointI * (views_.nTimeSteps() * fieldType_ + 3)
            + (timeI * fieldType_ + 3);
        for(label componentI = 0; componentI < fieldType_; componentI++)
        {
            reinterpret_cast<T1*>(values)[idx + componentI]
                =  component(meshMotion[pointI], componentI);
        }
    }
}

template <class T1>
void gmshViewPointMeshMotion<T1>
::getViewData(gmshViewBase::gmshViewBuffer& vB) const
{
    pointField meshMotion(mesh_.points().size(), vector::zero);

    gmshViewPointMesh<T1, vector>::
        writePointPoints(vB.buf[gmshViewBase::VP]);

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

        writePointMotions(vB.buf[gmshViewBase::VP], meshMotion, timeI);
    }
}

// * * * * * * * * * * * * * * * * Instantiators * * * * * * * * * * * * * * //

// double precision Gmsh view / pointMesh
template class gmshViewPointMeshMotion<double>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
