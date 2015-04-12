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
    Converts polyMesh to Gmsh mesh structure.

\*---------------------------------------------------------------------------*/

#include <fstream>
#include <iostream>

using std::ofstream;

#include "gmshElements.H"
#include "gmshMessageStream.H"
#include "polyMeshToGmsh.H"
#include "processorPolyPatch.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// physical region number starts from 1
label polyMeshToGmsh::totalPhysicals_ = 1;

    // Constructor
polyMeshToGmsh::polyMeshToGmsh(const IOobject& io, const label verbosity)
    : polyMesh(io), verbosity_(verbosity), nZoneFaces_(0), faceZoneI_(0),
      zoneFaceI_(0), physicalNameI_(0), vertexI_(0), elementI_(0),
      myProcNo_(0), myProcNoStr_()
{
    forAll(faceZones(), zoneI)
    {
        nZoneFaces_ += faceZones()[zoneI].size();
    }

    forAll(boundaryMesh(), patchI)
    {
        if(boundaryMesh().types()[patchI] == processorPolyPatch::typeName)
        {
            myProcNo_ = refCast<const processorPolyPatch>
                (boundaryMesh()[patchI]).myProcNo();
            OStringStream os;
            os << "(" << myProcNo_ << ")" ;
            myProcNoStr_ = os.str();
            break;
        }
    }
}

label polyMeshToGmsh::getNNames()
{
    // +1 for cells which does not belong to any zones
    return 1 + cellZones().size() + faceZones().size() + boundaryMesh().size();
}

label polyMeshToGmsh::getNElems()
{
    return nCells() + nZoneFaces_ + nFaces() - boundaryMesh()[0].start();
}

void polyMeshToGmsh::getPhysicalName(int& num, char name[])
{
    num = totalPhysicals_ + physicalNameI_;

    if(physicalNameI_ == 0) // for cells which does not belong to any zones
    {
        strncpy(name, (string("unnamed") + myProcNoStr_).c_str(), 256);
    }
    else if(physicalNameI_ <= cellZones().size())
    {
        strncpy(name,
        (cellZones().names()[physicalNameI_ - 1] + myProcNoStr_).c_str(), 256);
    }
    else if(physicalNameI_ <= cellZones().size() + faceZones().size())
    {
        const label zoneI = physicalNameI_ - 1 - cellZones().size();
        strncpy(name, (faceZones().names()[zoneI] + myProcNoStr_).c_str(),
        256);
    }
    else
    {
        const label patchI = physicalNameI_ - 1 - cellZones().size()
            - faceZones().size();
        std::string physicalName = boundaryMesh().names()[patchI];

        // omot processor number suffix for processor patches
        if(boundaryMesh().types()[patchI] != processorPolyPatch::typeName)
        {
            physicalName += myProcNoStr_;
        }

        // if the type is not patch append the type so that it can be retrieved
        // back by gmsh2ToFoam
        if(boundaryMesh().types()[patchI] != polyPatch::typeName)
        {
            physicalName += " " + boundaryMesh().types()[patchI];
        }

        strncpy(name, physicalName.c_str(), 256);
    }
    name[255] = '\0';

    physicalNameI_++;
}

void polyMeshToGmsh::getVertex(int& num, double& x, double& y, double& z)
{
    num = vertexI_;
    const pointField& p = points();
    x = p[vertexI_].x();
    y = p[vertexI_].y();
    z = p[vertexI_].z();
    vertexI_++;
}

void polyMeshToGmsh::getElementAttributes(int& num, int& type, int& physical,
int& elementary, int& partition) const
{
    static const cellModel& tet_ = *(cellModeller::lookup("tet"));
    static const cellModel& hex_ = *(cellModeller::lookup("hex"));
    static const cellModel& prism_ = *(cellModeller::lookup("prism"));
    static const cellModel& pyr_ = *(cellModeller::lookup("pyr"));

    if(elementI_ < nCells())
    {
        const cellModel& model = cellShapes()[elementI_].model();

        if(model == tet_)
        {
            type = 4;
        }
        else if(model == hex_)
        {
            type = 5;
        }
        else if(model == prism_)
        {
            type = 6;
        }
        else if(model == pyr_)
        {
            type = 7;
        }
        else
        {
            type = 0;
            gWarning << "Unhandled element type \"" << model.name()
                << "\" for element " << elementI_ << ": nPoints = "
                << model.nPoints() << ", nEdges = " << model.nEdges()
                << ", nFaces = " << model.nFaces() << endl;
        }

        num = elementI_;
        // whichZone returns -1 if the cell does not belong to any zones
        elementary = cellZones().whichZone(elementI_) + 1 + totalPhysicals_;
        physical = elementary;
        partition = myProcNo_;
    }
    else
    {
        label faceI;
        if(elementI_ < nCells() + nZoneFaces_)
        {
            faceI = faceZones()[faceZoneI_][zoneFaceI_];
            elementary = cellZones().size() + 1 + faceZoneI_ + totalPhysicals_;
        }
        else
        {
            faceI = elementI_ - nCells() - nZoneFaces_
                + boundaryMesh()[0].start();
            elementary = cellZones().size() + 1 + faceZones().size()
                + boundaryMesh().whichPatch(faceI) + totalPhysicals_;
        }

        const label nVerts = allFaces()[faceI].size();
        if(nVerts == 3)
        {
            type = 2;
        }
        else if(nVerts == 4)
        {
            type = 3;
        }
        else
        {
            type = 0;
            gWarning << "Unhandled number of face vertices "
                << nVerts << " for face " << faceI << endl;
        }

        num = elementI_;
        physical = elementary;
        partition = myProcNo_;
    }
}

void polyMeshToGmsh::getElementVerticesIndices(int indices[])
{
    if(elementI_ < nCells())
    {
        const labelList& cs = cellShapes()[elementI_];

        forAll(cs, pointI)
        {
            // Gmsh and polyMesh have the identical vertex numbering
            indices[pointI] = cs[pointI];
        }
    }
    else
    {
        label faceI;
        if(elementI_ < nCells() + nZoneFaces_)
        {
            faceI = faceZones()[faceZoneI_][zoneFaceI_];

            zoneFaceI_++;
            if(zoneFaceI_ >= faceZones()[faceZoneI_].size())
            {
                zoneFaceI_ = 0;
                faceZoneI_++;
            }
        }
        else
        {
            faceI = elementI_ - nCells() - nZoneFaces_
                + boundaryMesh()[0].start();
        }
        const face& f = allFaces()[faceI];

        forAll(f, pointI)
        {
            indices[pointI] = f[pointI];
        }
    }

    elementI_++;
}

void polyMeshToGmsh::writeGmshMesh()
{
    // make a directory called gmshInterface in the case
    mkDir(time().rootPath()/time().caseName()/"gmshInterface");

    fileName baseName(time().caseName());
    if(time().caseName() == ".")
    {
        baseName = time().rootPath().name();
    }
    else if(time().caseName() == "..")
    {
        baseName = time().rootPath().path().name();
    }

    // open a file for the mesh
    ofstream gmshMeshFile
    (
        (
            time().rootPath()/
            time().caseName()/
            "gmshInterface"/
            baseName + ".msh"
        ).c_str()
    );

    gInfo(verbosity_ >= 3) << "Writing Header" << endl;

    // write mesh format information
    gmshMeshFile << "$MeshFormat" << std::endl
        << "2.0 0 8" << std::endl
        << "$EndMeshFormat" << std::endl;

    // write physical names
    const label nNames = getNNames();
    if(nNames)
    {
        gmshMeshFile << "$PhysicalNames" << std::endl
            << nNames << std::endl;

        for(label nameI = 0; nameI < nNames; nameI++)
        {
            int num;
            char name[256];
            getPhysicalName(num, name);
            gmshMeshFile << num << " " << "\"" << name << "\"" << std::endl;
        }
        gmshMeshFile << "$EndPhysicalNames" << std::endl;
    }

    // write nodes
    const label nVerts = getNVerts();
    gmshMeshFile << "$Nodes" << std::endl
        << nVerts << std::endl;
    for(label vertI = 0; vertI < nVerts; vertI++)
    {
        int num;
        double xyz[3];
        getVertex(num, xyz[0], xyz[1], xyz[2]);
        gmshMeshFile << num << " " << xyz[0] << " " << xyz[1] << " " << xyz[2]
            << std::endl;
    }
    gmshMeshFile << "$EndNodes" << std::endl;

    // write elements
    const label nElems = getNElems();
    gmshMeshFile << "$Elements" << std::endl
        << nElems << std::endl;
    for(label elemI = 0; elemI < nElems; elemI++)
    {
        int num, type, physical, elementary, partition;
        getElementAttributes(num, type, physical, elementary, partition);
        gmshMeshFile << num << " " << type << " 3 " << physical << " "
            << elementary << " " << partition << " ";

        int indices[30];
        getElementVerticesIndices(indices);
        const label numVertices = gmshElements::typeToNVerts(type);
        for(label vertI = 0; vertI < numVertices - 1; vertI++)
        {
            gmshMeshFile << indices[vertI] << " ";
        }
        gmshMeshFile << indices[numVertices - 1] << std::endl;
    }
    gmshMeshFile << "$EndElements" << std::endl;
    gmshMeshFile.close();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
