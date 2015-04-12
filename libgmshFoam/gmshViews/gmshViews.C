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
    Reads OpenFOAM fields and converts to the Gmsh View format.

\*---------------------------------------------------------------------------*/

#include "IOobjectList.H"
#include "cellModeller.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"

#include "gmshViews.H"
#include "gmshView.H"
#include "gmshViewMeshMotion.H"
#include "gmshMessageStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Static data   * * * * * * * * * * * * * * //

static const char *fieldClasses[]
= {"volScalarField", "volVectorField", "volTensorField",
   "surfaceScalarField", "surfaceVectorField", "surfaceTensorField"};

const wordListFromCharArray gmshViews::fieldClasses_(6, fieldClasses);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gmshViews::gmshViews(bool& good, const fileName& rootPath,
const fileName& caseName, const gmshViewsOptions& opt)
    :
    runTime_(Foam::Time::controlDictName, rootPath, caseName),
    // construct mesh
    mesh_
    (
        IOobject
        (
            fvMesh::defaultRegion, runTime_.timeName(),
            runTime_, IOobject::MUST_READ
        )
    ),
    timeList_(runTime_.times()),
    fieldNames_(0),
    // get the startTime_, and also initializes fieldNames_
    startTime_(findStartTime(good, opt)),
    nTimes_(timeList_.size()),
    fieldI_(0),
    meshTimeI_(0),
    isMeshMotion_(false),
    isVolMeshMotionDone_(false),
    isSurfaceMeshMotionDone_(false),
    motionClassI_(0),
    gV_(NULL),
    splitTimeStepsByMeshMotion_(opt.splitTimeStepsByMeshMotion_),
    verbosity_(opt.verbosity_)
{
    // mesh.constructAndClear();

    if(!good)
    {
        gSeriousError << "Could not find supported field objects." << endl;
        return;
    }

    // print the fields found
    gInfo(verbosity_ >= 3) << endl << "Found supported objects at time = "
        << timeList_[startTime_].name() << ":" << endl;
    forAll(fieldNames_, classI)
    {
        if(fieldNames_[classI].size() || verbosity_ >= 4)
        {
            gInfo(verbosity_ >= 3) << "    " << fieldClasses_[classI] << ": ";
            forAll(fieldNames_[classI], fieldI)
            {
                gInfo(verbosity_ >= 3) << fieldNames_[classI][fieldI] << " ";
            }
            gInfo(verbosity_ >= 3) << endl;
        }
    }

    gInfo(splitTimeStepsByMeshMotion_ && verbosity_ >= 3) << endl
        << "Splitting field views by mesh motion:" << endl;

    // detect mesh motions
    isMeshMotion_ = false;
    meshTime_.append(startTime_);
    for(label timeI = startTime_; timeI < timeList_.size(); timeI++)
    {
        runTime_.setTime(timeList_[timeI], timeI);
        IOobject ioPoints("points", runTime_.timeName(), polyMesh::typeName,
        mesh_, IOobject::MUST_READ, IOobject::NO_WRITE);

        if ((runTime_.timeName() != runTime_.constant())
        && ioPoints.headerOk())
        {
            if(splitTimeStepsByMeshMotion_)
            {
                meshTime_.append(timeI);
                gInfo(verbosity_ >= 3) << "    Found mesh motion at t = "
                    << timeList_[timeI].name() << endl;
            }
            else
            {
                isMeshMotion_ = true;
                gInfo(verbosity_ >= 3) << endl << "Detected mesh motion at t = "
                    << timeList_[timeI].name()
                    << ". Constructing mesh motion fields." << endl;
                break;
            }
        }
    }
    meshTime_.append(timeList_.size());
    meshTime_.shrink();
    endTime_ = meshTime_[1] - 1;
    nTimeSteps_ = endTime_ - startTime_ + 1;

    gInfo(verbosity_ >= 3) << "    nTimeSteps = " << nTimeSteps_ << endl;

    // construct the element lists
    makeFaceLists();
    makeShapeCellLists();
    makeBoundaryFaceLists();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// find the startTime_ and sets fieldNames_
label gmshViews::findStartTime(bool& good, const gmshViewsOptions& opt)
{
    label startTime = 0;
    good = false;

    fieldNames_.setSize(fieldClasses_.size());

    // determine the starting time
    label timeI;
    if(opt.startTimeStr_ == "")
    {
        timeI = 0;
    }
    else if(opt.startTimeStr_ == "latestTime")
    {
        timeI = timeList_.size() - 1;
    }
    else
    {
        for(timeI = 0; timeI < timeList_.size(); timeI++)
        {
            if(opt.startTimeStr_ == timeList_[timeI].name())
            {
                break;
            }
        }
        if(timeI == timeList_.size())
        {
            gSeriousError(opt.verbosity_ >= 1)
                << "Could not find the time with name \""
                    << opt.startTimeStr_.c_str() << "\"." << endl;
            return timeI; // with good == false
        }
    }

    gInfo(opt.verbosity_ >= 4)
        << "Starting searching field objects from time = "
            << timeList_[timeI].name() << endl;

    // search at a time for the given field classes
    for(; timeI < timeList_.size(); timeI++)
    {
        runTime_.setTime(timeList_[timeI], timeI);
        IOobjectList objects(mesh_, runTime_.timeName());

        // for all the supported field classes:
        forAll(fieldClasses_, classI)
        {
            // if a field is found at a time the starting times for all the
            // fields are fixed to the time.
            if(objects.names(fieldClasses_[classI]).size())
            {
                good = true;
                startTime = timeI;
                classI_ = classI;
                startClass_ = classI;

                // copy all the field names which belong to all the
                // remaining classes
                for(label classJ = classI; classJ < fieldClasses_.size();
                    classJ++)
                {
                    // count the number of fields for the class
                    label nFields
                        = objects.names(fieldClasses_[classJ]).size();
                    // copy field names for classJ
                    fieldNames_[classJ].setSize(nFields);
                    forAll(fieldNames_[classJ], fieldI)
                    {
                        fieldNames_[classJ][fieldI]
                            = objects.names(fieldClasses_[classJ])[fieldI];
                    }
                }
                break;
            }
            else
            {
                fieldNames_[classI].setSize(0);
            }
        }
        if(good)
        {
            break;
        }
    }
    return startTime;
}

void gmshViews::makeBoundaryFaceLists()
{
    label nPatches = mesh_.boundaryMesh().size();
    patchTris_.setSize(nPatches);
    patchQuads_.setSize(nPatches);

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];

        if(isType<emptyPolyPatch>(pp))
        {
            gInfo(verbosity_ >= 4)
                << "    Skipping the " << mesh_.boundaryMesh().types()[patchI]
                    << " patch " << mesh_.boundaryMesh().names()[patchI] << "."
                    << endl;
        }
        else
        {
            forAll(pp, faceI)
            {
                const face& f = pp[faceI];
                if(f.size() == 3)
                {
                    patchTris_[patchI].append(faceI);
                }
                else if(f.size() == 4)
                {
                    patchQuads_[patchI].append(faceI);
                }
                else
                {
                    gWarning << "Unhandled number of face vertices "
                        << f.size() << " for face " << faceI << " of patch "
                        << patchI << endl;
                }
            }
        }
        patchTris_[patchI].shrink();
        patchQuads_[patchI].shrink();
    }
}

void gmshViews::makeFaceLists()
{
    for(label faceI = 0; faceI < mesh_.boundaryMesh()[0].start(); faceI++)
    {
        const face& f = mesh_.faces()[faceI];

        if(f.size() == 3)
        {
            tris_.append(faceI);
        }
        else if(f.size() == 4)
        {
            quads_.append(faceI);
        }
        else
        {
            gWarning << "Unhandled number of face vertices " << f.size()
                << " for face " << faceI << endl;
        }
    }
    tris_.shrink();
    quads_.shrink();
}

void gmshViews::makeShapeCellLists()
{
    const cellModel& tet = *(cellModeller::lookup("tet"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& hex = *(cellModeller::lookup("hex"));

    const cellShapeList& cellShapes = mesh_.cellShapes();

    forAll(cellShapes, cellI)
    {
        const cellModel& model = cellShapes[cellI].model();

        if(model == tet)
        {
            tets_.append(cellI);
        }
        else if(model == hex)
        {
            hexes_.append(cellI);
        }
        else if(model == prism)
        {
            prisms_.append(cellI);
        }
        else if(model == pyr)
        {
            pyrs_.append(cellI);
        }
        else
        {
            gWarning << "Unhandled element type \"" << model.name()
                << "\" for cell " << cellI << ": nPoints = " << model.nPoints()
                << ", nEdges = " << model.nEdges() << ", nFaces = "
                << model.nFaces() << endl;
        }
    }
    tets_.shrink();
    hexes_.shrink();
    prisms_.shrink();
    pyrs_.shrink();
}

gmshViewBase *gmshViews::getNextView()
{
    if(gV_ != NULL)
    {
        delete gV_;

        fieldI_++;
        if(fieldI_ >= fieldNames_[classI_].size())
        {
            fieldI_ = 0;

            // skip the field classes where no corresponding field is present
            for(classI_++;
                classI_ < fieldNames_.size()
                    && fieldNames_[classI_].size() == 0;
                classI_++);

            if(classI_ >= fieldNames_.size())
            {
                if(splitTimeStepsByMeshMotion_)
                {
                    meshTimeI_++;
                    if(meshTimeI_ + 1 >= meshTime_.size())
                    {
                        gV_ = NULL;
                        return gV_;
                    }

                    startTime_ = meshTime_[meshTimeI_];
                    endTime_ = meshTime_[meshTimeI_ + 1] - 1;
                    nTimeSteps_ = endTime_ - startTime_ + 1;

                    runTime_.setTime(timeList_[startTime_], startTime_);
                    IOobject ioPoints("points", runTime_.timeName(),
                    polyMesh::typeName, mesh_);

                    if(ioPoints.headerOk())
                    {
                        gInfo(verbosity_ >= 3) << endl
                            << "Applying mesh movement at t = "
                            << timeList_[startTime_].name() << endl;

                        pointIOField newPoints(ioPoints);
                        mesh_.movePoints(newPoints);
                    }
                    else
                    {
                        gSeriousError << "Failed reading new mesh points at "
                            << "t = " << runTime_.timeName() << endl;
                    }

                    classI_ = startClass_;
                    fieldI_ = 0;
                }
                else
                {
                    gV_ = NULL;
                    return gV_;
                }
            }
        }
    }

    if(fieldClasses_[classI_] == "volScalarField")
    {
        gV_ = new gmshView<double, scalar, volMesh>(*this,
        fieldNames_[classI_][fieldI_], runTime_, mesh_, verbosity_);
    }
    else if(fieldClasses_[classI_] == "volVectorField")
    {
        gV_ = new gmshView<double, vector, volMesh>(*this,
        fieldNames_[classI_][fieldI_], runTime_, mesh_, verbosity_);
    }
    else if(fieldClasses_[classI_] == "volTensorField")
    {
        gV_ = new gmshView<double, tensor, volMesh>(*this,
        fieldNames_[classI_][fieldI_], runTime_, mesh_, verbosity_);
    }
    else if(fieldClasses_[classI_] == "surfaceScalarField")
    {
        gV_ = new gmshView<double, scalar, surfaceMesh>(*this,
        fieldNames_[classI_][fieldI_], runTime_, mesh_, verbosity_);
    }
    else if(fieldClasses_[classI_] == "surfaceVectorField")
    {
        gV_ = new gmshView<double, vector, surfaceMesh>(*this,
        fieldNames_[classI_][fieldI_], runTime_, mesh_, verbosity_);
    }
    else if(fieldClasses_[classI_] == "surfaceTensorField")
    {
        gV_ = new gmshView<double, tensor, surfaceMesh>(*this,
        fieldNames_[classI_][fieldI_], runTime_, mesh_, verbosity_);
    }
    else
    {
        gSeriousError << "Unsupported field class " << fieldClasses_[classI_]
            << endl;
        gV_ = NULL;
        return NULL;
    }

    return gV_;
}

gmshViewBase *gmshViews::getNextMeshMotion()
{
    if(gV_ != NULL)
    {
        delete gV_;
        gV_ = NULL;
    }

    if(!isMeshMotion_)
    {
        return NULL;
    }

    while(motionClassI_ < fieldClasses_.size())
    {
        if(fieldClasses_[motionClassI_].substr(0, 3) == "vol"
        && fieldNames_[motionClassI_].size() && !isVolMeshMotionDone_)
        {
            gV_ = new gmshViewMeshMotion<double, volMesh>(*this,
            "volMeshMotion", runTime_, mesh_, verbosity_);
            isVolMeshMotionDone_ = true;
            motionClassI_++;
            return gV_;
        }
        else if(fieldClasses_[motionClassI_].substr(0, 7) == "surface"
        && fieldNames_[motionClassI_].size() && !isSurfaceMeshMotionDone_)
        {
            gV_ = new gmshViewMeshMotion<double, surfaceMesh>(*this,
            "surfaceMeshMotion", runTime_, mesh_, verbosity_);
            isSurfaceMeshMotionDone_ = true;
            motionClassI_++;
            return gV_;
        }
        motionClassI_++;
    }

    return NULL;
}

void gmshViews::convertView(gmshViewBase& gV, const fileName& gmshViewsFolder)
{
    // get the format
    gmshViewBase::gmshPostFormat pF;
    gV.getPostFormat(pF);

    // get the header information
    gmshViewBase::gmshViewHeader vH;
    gV.getViewHeader(vH);

    // allocate memory
    gmshViewBase::gmshViewBuffer vB;
    label sSize[gmshViewBase::nGeoTypes];
    gV.getViewSize(vH, sSize);
    for(label typeI = 0; typeI < gmshViewBase::nGeoTypes; typeI++)
    {
        vB.buf[typeI]
            = new char[sSize[typeI] ? pF.size * sSize[typeI] : 1];
    }
    vB.T2D = new char[vH.NbT2 ? 4 * vH.NbT2 * pF.size : 1];
    vB.T2C = new char[vH.t2l ? vH.t2l * sizeof(char) : 1];
    vB.T3D = new char[vH.NbT3 ? 5 * vH.NbT3 * pF.size : 1];
    vB.T3C = new char[vH.t3l ? vH.t3l * sizeof(char) : 1];

    // get timestep values
    char *tB = new char [pF.size * vH.NbTimeStep];
    gV.getTimeStepValues(tB);

    // get view data
    gV.getViewData(vB);

    // write view
    if(splitTimeStepsByMeshMotion_)
    {
        gV.write(gmshViewsFolder/vH.name + "_" + timeList_[startTime_].name()
        + ".pos", pF, vH, tB, vB);
    }
    else
    {
        gV.write(gmshViewsFolder/vH.name + ".pos", pF, vH, tB, vB);
    }

    // clear buffers
    delete [] tB;
    for(label typeI = 0; typeI < gmshViewBase::nGeoTypes; typeI++)
    {
        delete [] vB.buf[typeI];
    }
    delete [] vB.T2D;
    delete [] vB.T2C;
    delete [] vB.T3D;
    delete [] vB.T3C;
}

void gmshViews::convert(const fileName& gmshViewsFolder)
{
    gmshViewBase* gVPtr;

    // get the next mesh motion
    while((gVPtr = getNextMeshMotion()) != NULL)
    {
        convertView(*gVPtr, gmshViewsFolder);
    }

    // get the next view
    while((gVPtr = getNextView()) != NULL)
    {
        convertView(*gVPtr, gmshViewsFolder);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
