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
    gmshElements

Description
    The class holds the $Elements section of a Gmsh .msh file.

SourceFiles
    gmshElements.C

\*---------------------------------------------------------------------------*/

#ifndef gmshElements_H
#define gmshElements_H

#include "IFstream.H"
#include "cellModeller.H"
#include "polyMesh.H"

namespace Foam
{

// The $Elements section data structure
class gmshElements
{
    // Private structure
    struct elementTypes {
        bool isHandled_;
        label nVerts_;
        const cellModel* cellModel_;
        char *name_;
    };
    const static struct elementTypes elemTab_[];
    const static label maxNVerts;

    const label verbosity_;

    // Storage for all cells.
    cellShapeList cells_;

    // Region number <--> patch/zone mappings
    Map<label> regionToPatch_;
    Map<label> regionToZone_;
    DynamicList<label> patchToRegion_;
    DynamicList<label> zoneToRegion_;

    // Patch/zone lists
    List<DynamicList<face> > patchFaces_;
    List<DynamicList<label> > zoneCells_;

    // Class-global variables

    label cellI_;
    labelList nCells_;

    // Is physical region defined or not?
    // If $PhysicalNames section is present it is obviously defined
    bool isPhysicalRegion_;

    // Reference to Gmsh node numbers -> OpenFOAM point numbers mapper
    Map<label>& mshToFoam_;

    // From gmsh physical / elementary region to Foam patch
    Map<label> physicalRegionToPatch_, elementaryRegionToPatch_;
    // From gmsh physical / elementary region to Foam cellZone
    Map<label> physicalRegionToZone_, elementaryRegionToZone_;

    // Reverse mapping for (physical/elementary)RegionToPatch
    DynamicList<label> elementaryPatchToRegion_, physicalPatchToRegion_;
    // Reverse mapping for (physical/elementary)ZoneToRegion
    DynamicList<label> elementaryZoneToRegion_, physicalZoneToRegion_;

    // Storage for patch faces.
    List<DynamicList<face> > physicalPatchFaces_, elementaryPatchFaces_;
    // Storage for cell zones.
    List<DynamicList<label> > physicalZoneCells_, elementaryZoneCells_;

    // Private functions
    void renumberPoints(const Map<label>& mshToFoam, labelList& labels);

    void storeFaceInPatch(const label reg, const face& polygon,
    Map<label>& regionToPatch, List<DynamicList<face> >& patchFaces,
    DynamicList<label>& patchToRegion, const word &physicalOrElementary);

    void storeCellInZone(const label regPhys, const label cellI,
    Map<label>& regionToZone, List<DynamicList<label> >& zoneCells,
    DynamicList<label>& zoneToRegion, const word& physicalOrElementary);

public:

    gmshElements(Map<label>& mshToFoam, const label verbosity):
        verbosity_(verbosity), mshToFoam_(mshToFoam)
	{}
    const cellShapeList& cells() const { return cells_; }
    cellShapeList& cells() { return cells_; }
    const Map<label>& regionToPatch() const { return regionToPatch_; }
    Map<label>& regionToPatch() { return regionToPatch_; }
    const Map<label>& regionToZone() const { return regionToZone_; }
    Map<label>& regionToZone() { return regionToZone_; }
    const DynamicList<label>& patchToRegion() const { return patchToRegion_; }
    const DynamicList<label>& zoneToRegion() const { return zoneToRegion_; }
    const List<DynamicList<face> >& patchFaces() const { return patchFaces_; }
    List<DynamicList<face> >& patchFaces() { return patchFaces_; }
    const List<DynamicList<label> >& zoneCells() const { return zoneCells_; }

    void readElements(IFstream& inFile, const Map<label>& mshToFoam,
    const bool isPhysicalNames, const scalar formatVersion,
    const word& execName);

    // for Gmsh direct interface
    void setNElems(const label nElems, const bool isPhysicalRegion);
    void insertElement(const label elmNumber, const label elmType,
    const label regPhys, const label regElem, const label partition,
    const label vertices[]);
    void postInsertElements();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
