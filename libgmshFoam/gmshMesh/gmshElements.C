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
    Reads the $Elements section of a .msh file.

\*---------------------------------------------------------------------------*/

#include "gmshElements.H"
#include "gmshMessageStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Static data   * * * * * * * * * * * * * * //

const label gmshElements::maxNVerts = 27;
const struct gmshElements::elementTypes gmshElements::elemTab_[] = {
    {false, 0, NULL, ""},                                // 0:(dummy)
    {false, 2, NULL, "   line"},                         // 1:LINE
    {true, 3, NULL, "    tri"},                          // 2:TRI
    {true, 4, NULL, "   quad"},                          // 3:QUAD
    {true, 4, cellModeller::lookup("tet"), "    tet"},   // 4:TET
    {true, 8, cellModeller::lookup("hex"), "    hex"},   // 5:HEX
    {true, 6, cellModeller::lookup("prism"), "  prism"}, // 6:PRISM
    {true, 5, cellModeller::lookup("pyr"), "    pyr"},   // 7:PYR
    {false, 3, NULL, "  line2"},                         // 8:LINE2
    {false, 6, NULL, "   tri2"},                         // 9:TRI2
    {false, 9, NULL, "  quad2"},                         //10:QUAD2
    {false, 10, NULL, "   tet2"},                        //11:TET2
    {false, 27, NULL, "   hex2"},                        //12:HEX2
    {false, 18, NULL, " prism2"},                        //13:PRISM2
    {false, 14, NULL, "   pyr2"},                        //14:PYR2
    {false, 1, NULL, "  point"},                         //15:PNT
    {false, 8, NULL, " quad2i"},                         //16:QUAD2I
    {false, 20, NULL, "  hex2i"},                        //17:HEX2I
    {false, 15, NULL, "prism2i"},                        //18:PRISM2I
    {false, 13, NULL, "  pyr2i"}                         //19:PYR2I
};

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gmshElements::renumberPoints(const Map<label>& mshToFoam,
labelList& labels)
{
    forAll(labels, labelI)
    {
        labels[labelI] = mshToFoam[labels[labelI]];
    }
}

void gmshElements::storeFaceInPatch(const label reg, const face& polygon,
Map<label>& regionToPatch, List<DynamicList<face> >& patchFaces,
DynamicList<label>& patchToRegion, const word &physicalOrElementary)
{
    // Find the patch where the element belongs
    Map<label>::iterator regFnd = regionToPatch.find(reg);
    label patchI = -1;
    if (regFnd == regionToPatch.end())
    {
        // New region. Allocate patch for it.
        patchFaces.setSize(patchFaces.size() + 1);

        patchI = patchFaces.size()-1;

        gInfo(verbosity_ >= 5) << "    Mapping " << physicalOrElementary
            << " region " << reg << " to Foam " << physicalOrElementary
            << " patch/faceZone " << patchI << endl;

        regionToPatch.insert(reg, patchI);

        patchToRegion.append(reg);
    }
    else
    {
        // Existing patch for region
        patchI = regFnd();
    }

    // Add polygon to correct PatchFaces.
    patchFaces[patchI].append(polygon);
}

void gmshElements::storeCellInZone(const label regPhys, const label cellI,
Map<label>& regionToZone, List<DynamicList<label> >& zoneCells,
DynamicList<label>& zoneToRegion, const word& physicalOrElementary)
{
    Map<label>::iterator zoneFnd = regionToZone.find(regPhys);

    if (zoneFnd == regionToZone.end())
    {
        // New region. Allocate zone for it.
        zoneCells.setSize(zoneCells.size() + 1);

        label zoneI = zoneCells.size()-1;

        gInfo(verbosity_ >= 5) << "    Mapping " << physicalOrElementary
            << " region " << regPhys << " to Foam " << physicalOrElementary
            << " cellZone " << zoneI << endl;

        regionToZone.insert(regPhys, zoneI);

        zoneCells[zoneI].append(cellI);

	zoneToRegion.append(regPhys);
    }
    else
    {
        // Existing zone for region
        zoneCells[zoneFnd()].append(cellI);
    }
}

void gmshElements::setNElems(const label nElems, const bool isPhysicalRegion)
{
    // cells_ too big. Shrink later
    cells_.setSize(nElems);

    nCells_.setSize(20, 0);

    // initial guess of physical region definition
    isPhysicalRegion_ = isPhysicalRegion;

    cellI_ = 0;
}

void gmshElements::insertElement(const label elmNumber, const label elmType,
const label regPhys, const label regElem, const label partition,
const label vertices[])
{
    // partition information is not used for now

    if(elmType < 1 || elmType > 19)
    {
        gWarning << "Undefined element type " << elmType
            << " in element number " << elmNumber << endl;
        return;
    }

    // if a value other than zero is found in regPhys the physical region
    // is determined to be defined
    if(isPhysicalRegion_ == false && regPhys != 0)
    {
        isPhysicalRegion_ = true;

        Info
            << "Found a physical region of region number other than 0; "
            << endl << "omitting elementary region handling." << endl;

        elementaryRegionToPatch_.clear();
        elementaryRegionToZone_.clear();
        elementaryPatchToRegion_.clear();
        elementaryZoneToRegion_.clear();
        elementaryPatchFaces_.clear();
        elementaryZoneCells_.clear();
    }

    if(elemTab_[elmType].isHandled_)
    {
        if(elemTab_[elmType].cellModel_ == NULL)
        {
            // The element is a face
            face facePoints(elemTab_[elmType].nVerts_);

            forAll(facePoints, vertI)
            {
                facePoints[vertI] = vertices[vertI];
            }
            renumberPoints(mshToFoam_, facePoints);

            // Find the physical patch where the element belongs
            storeFaceInPatch(regPhys, facePoints, physicalRegionToPatch_,
            physicalPatchFaces_, physicalPatchToRegion_, "physical");

            if(isPhysicalRegion_ == false)
            {
                // Do the same thing for the elementary patch
                storeFaceInPatch(regElem, facePoints,
                elementaryRegionToPatch_, elementaryPatchFaces_,
                elementaryPatchToRegion_, "elementary");
            }

            nCells_[elmType]++;
        }
        else
        {
            // The element is a cell
            storeCellInZone(regPhys, cellI_, physicalRegionToZone_,
            physicalZoneCells_, physicalZoneToRegion_, "physical");

            if(isPhysicalRegion_ == false)
            {
                storeCellInZone(regElem, cellI_, elementaryRegionToZone_,
                elementaryZoneCells_, elementaryZoneToRegion_, "elementary");
            }

            labelList cellPoints(elemTab_[elmType].nVerts_);

            forAll(cellPoints, vertI)
                {
                    cellPoints[vertI] = vertices[vertI];
                }

            renumberPoints(mshToFoam_, cellPoints);

            cells_[cellI_++]
                = cellShape(*elemTab_[elmType].cellModel_, cellPoints);

            nCells_[elmType]++;
        }
    }
    else
    {
        nCells_[elmType]++;
    }
}

void gmshElements::postInsertElements()
{
    cells_.setSize(cellI_);

    if(isPhysicalRegion_)
    {
        regionToPatch_ = physicalRegionToPatch_;
        regionToZone_ = physicalRegionToZone_;
        patchFaces_ = physicalPatchFaces_;
        zoneCells_ = physicalZoneCells_;
        patchToRegion_ = physicalPatchToRegion_;
        zoneToRegion_ = physicalZoneToRegion_;
    }
    else
    {
        regionToPatch_ = elementaryRegionToPatch_;
        regionToZone_ = elementaryRegionToZone_;
        patchFaces_ = elementaryPatchFaces_;
        zoneCells_ = elementaryZoneCells_;
        patchToRegion_ = elementaryPatchToRegion_;
        zoneToRegion_ = elementaryZoneToRegion_;

        elementaryRegionToPatch_.clear();
        elementaryRegionToZone_.clear();
        elementaryPatchFaces_.clear();
        elementaryZoneCells_.clear();
        elementaryPatchToRegion_.clear();
        elementaryZoneToRegion_.clear();
    }

    physicalRegionToPatch_.clear();
    physicalRegionToZone_.clear();
    physicalPatchFaces_.clear();
    physicalZoneCells_.clear();
    physicalPatchToRegion_.clear();
    physicalZoneToRegion_.clear();

    forAll(patchFaces_, patchI)
    {
        patchFaces_[patchI].shrink();
    }

    forAll(zoneCells_, zoneI)
    {
        zoneCells_[zoneI].shrink();
    }

    gInfo(verbosity_ >= 4) << endl << "Cells:" << endl << "      total: "
        << cells_.size() << endl;
    for(label typeI = 1; typeI <= 19; typeI++)
    {
        if(verbosity_ >= 5 || (verbosity_ >= 4 && nCells_[typeI] > 0))
        {
            Info << "    " << elemTab_[typeI].name_ << ": " << nCells_[typeI];
            if(!elemTab_[typeI].isHandled_)
            {
                Info << " (not converted to OpenFOAM mesh)";
            }
            else if(elemTab_[typeI].cellModel_ == NULL)
            {
                Info << " (not counted as cells)";
            }
            Info << endl;
        }
    }
}


void gmshElements::readElements(IFstream& inFile, const Map<label>& mshToFoam,
const bool isPhysicalNames, const scalar formatVersion, const word& execName)
{
    label nElems;
    inFile >> nElems;

    Info << endl << "Read nElems: " << nElems << endl;

    setNElems(nElems, isPhysicalNames);

    gInfo(verbosity_ >= 4) << endl
        << "Reading elements and mapping Gmsh regions to Foam patches/zones"
        << (verbosity_ >= 5 ? ":" : "") << endl;

    forAll(cells_, elemI)
    {
        gInfo(verbosity_ >= 5 && (elemI % 1000) == 0)
            << "        Reading element " << elemI << " out of "
                << cells_.size() << endl;

        label elmNumber, elmType, regPhys, regElem, partition = 0;

	if(formatVersion == 1.0)
        {
	    label nNodes;
	    inFile >> elmNumber >> elmType >> regPhys >> regElem >> nNodes;
        }
	else if(formatVersion == 2.0)
        {
	    label nTags;
	    inFile >> elmNumber >> elmType >> nTags;
	    if (nTags < 2)
            {
		FatalErrorIn(execName)
                    << "Too few tags on line " << inFile.lineNumber()
                        << exit(FatalError);
            }
	    inFile >> regPhys >> regElem;

            label tagI;
            if(nTags >= 3)
            {
                inFile >>  partition;
                tagI = 3;
            }
            else
            {
                partition = 0;
                tagI = 2;
            }
	    for(; tagI < nTags; tagI++)
            {
		label tagInt;
		inFile >> tagInt;
            }
        }

        label vertices[maxNVerts];
        for(label vertI = 0; vertI < elemTab_[elmType].nVerts_; vertI++)
        {
            inFile >> vertices[vertI];
        }

        insertElement(elmNumber, elmType, regPhys, regElem, partition,
        vertices);
    }

    word tag(inFile);
    if ((formatVersion == 1.0 && tag != "$ENDELM" && tag != "$ENDELM\r")
    || (formatVersion == 2.0 && tag != "$EndElements"
    && tag != "$EndElements\r"))
    {
        FatalErrorIn(execName)
            << "Did not find $ENDNOD nor $EndNodes tags on line "
                << inFile.lineNumber() << exit(FatalError);
    }

    gInfo(!isPhysicalRegion_ && verbosity_ >= 3)
        << "Physical region of region number other than 0 not found;"
            << endl << "using elementary region(s) instead." << endl;

    postInsertElements();

    // All OpenFOAM meshes must be three-dimensional.
    if(cells_.size() == 0)
    {
        FatalErrorIn(execName)
            << "No volumetric cells (tets/pyrs/prisms/hexes) found" << nl
                << "    in " << inFile.name() << ". Aborting conversion."
                << " Check your mesh carefully," << nl
                << "    particularly whether your physical volume definition"
                << " is correct."
                << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
