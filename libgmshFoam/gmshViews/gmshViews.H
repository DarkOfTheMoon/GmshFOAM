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
    gmshViews

Description
    Reads OpenFOAM fields and converts to the Gmsh View format.

SourceFiles
    gmshViews.C

\*---------------------------------------------------------------------------*/

#ifndef gmshViews_H
#define gmshViews_H

#ifdef cygwin
#include "Time.hh"
#else
#include "Time.H"
#endif
#include "fvMesh.H"

#include "gmshViewBase.H"
#include "wordListFromCharArray.H"

namespace Foam
{

class gmshViews
{
public:

    struct gmshViewsOptions
    {
        label nOptions_;
        label splitTimeStepsByMeshMotion_;
        label verbosity_;
        string startTimeStr_;
    };

private:

    static const wordListFromCharArray fieldClasses_;

    Time runTime_;
    fvMesh mesh_;
    const instantList timeList_;
    List<wordList> fieldNames_;
    label startTime_, endTime_;
    const label nTimes_;
    label nTimeSteps_;

    label startClass_, classI_;
    label fieldI_;
    label meshTimeI_;
    bool isMeshMotion_;
    bool isVolMeshMotionDone_, isSurfaceMeshMotionDone_,
        isPointMeshMotionDone_;
    label  motionClassI_;

    gmshViewBase *gV_;

    // Mesh properties
    DynamicList<label> tris_, quads_, tets_, hexes_, prisms_, pyrs_;
    List<DynamicList<label> > patchTris_, patchQuads_;
    DynamicList<label> meshTime_;

    // Split timesteps of a field view by mesh motions
    const bool splitTimeStepsByMeshMotion_;
    // Verbosity
    const label verbosity_;

    // Private member functions
    label findStartTime(bool& good, const gmshViewsOptions& opt);
    void makeFaceLists();
    void makeBoundaryFaceLists();
    void makeShapeCellLists();

public:

    // Constructor
    gmshViews(bool& good, const fileName& rootPath, const fileName& caseName,
    const gmshViewsOptions& opt);

    const instantList& timeList() const { return timeList_; }
    const label startTime() const { return startTime_; }
    const label endTime() const { return endTime_; }
    const label nTimes() const { return nTimes_; }
    const label nTimeSteps() const { return nTimeSteps_; }
    const bool isMeshMotion() const { return isMeshMotion_; }

    const DynamicList<label>& tris() const { return tris_; }
    const DynamicList<label>& quads() const { return quads_; }
    const labelList& tets() const { return tets_; }
    const labelList& hexes() const { return hexes_; }
    const labelList& prisms() const { return prisms_; }
    const labelList& pyrs() const { return pyrs_; }

    const List<DynamicList<label> >& patchTris() const { return patchTris_; }
    const List<DynamicList<label> >& patchQuads() const { return patchQuads_; }

    gmshViewBase *getNextMeshMotion();
    gmshViewBase *getNextView();
    void convertView(gmshViewBase& gV, const fileName& gmshViewsFolder);
    void convert(const fileName& gmshViewsFolder);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
