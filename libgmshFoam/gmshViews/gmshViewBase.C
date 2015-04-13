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
    The abstract base class for reading an OpenFOAM field and converts
    to Gmsh Views.

\*---------------------------------------------------------------------------*/

#include "OFstream.H"

#include "gmshMessageStream.H"
#include "gmshViewBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Static data * * * * * * * * * * * * * * * //

label gmshViewBase::nVertices_[]
    = {1, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14};
label gmshViewBase::nComponents_[] = {1, 3, 9};

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gmshViewBase::getViewSize(const gmshViewHeader& vH, label vSize[]) const
{
    for(label typeJ = 0; typeJ < 15; typeJ++)
    {
        for(label typeI = 0; typeI < 3; typeI++)
        {
            vSize[typeJ * 3 + typeI] =
                vH.Nb[typeJ * 3 + typeI] *
                (
                    nVertices_[typeJ] *
                    (
                        vH.NbTimeStep * nComponents_[typeI] + 3
                    )
                );
        }
    }
}

void gmshViewBase::write(const fileName& name, const gmshPostFormat& pF,
const gmshViewHeader& vH, const char* tB, const gmshViewBuffer &vB) const
{
    // using std::ofstream because Foam::OFstream writes unnecessary delimiters
    std::ofstream vFile(name.c_str(), std::ios::binary);
    OStringStream oStr;

    Info << "    Writing view" << endl;
    gInfo(verbosity_ >= 4) << "        Writing header" << endl;

    // write header
    oStr << "$PostFormat" << endl
        << pF.version << " " << pF.format << " " <<  pF.size << endl
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

    if(pF.format == 1)
    {
        // format is binary

        const int one = 1;
        vFile.write(reinterpret_cast<const char*>(&one), sizeof(int));

        vFile.write(tB, pF.size * vH.NbTimeStep);

        gInfo(verbosity_ >= 4) << "        Writing elements" << endl;

        label sSize[gmshViewBase::nGeoTypes];
        getViewSize(vH, sSize);

        for(label typeI = 0; typeI < gmshViewBase::nGeoTypes; typeI++)
        {
            vFile.write(vB.buf[typeI], pF.size * sSize[typeI]);
        }

        // 2D texts
        vFile.write(vB.T2D, vH.NbT2 * 4 * pF.size);
        vFile.write(vB.T2C, vH.t2l * sizeof(char));

        // 3D texts
        vFile.write(vB.T3D, vH.NbT3 * 5 * pF.size);
        vFile.write(vB.T3C, vH.t3l * sizeof(char));
    }
    else
    {
        // format is ascii

        for(label timeI = 0; timeI < vH.NbTimeStep; timeI++)
        {
            if(pF.size == 4)
            {
                vFile << reinterpret_cast<const float*>(tB)[timeI] << " ";
            }
            else
            {
                vFile << reinterpret_cast<const double*>(tB)[timeI] << " ";
            }
        }

        gInfo(verbosity_ >= 4) << "        Writing elements" << endl;

        vFile << std::endl;
        for(label typeI = 0; typeI < gmshViewBase::nGeoTypes; typeI++)
        {
            label elemI = 0;
            for(label cellI = 0; cellI < vH.Nb[typeI]; cellI++)
            {
                for(label componentI = 0; componentI < 3; componentI++)
                {
                    for(label pointI = 0; pointI < nVertices_[typeI / 3];
                        pointI++)
                    {
                        if(pF.size == 4)
                        {
                            vFile << reinterpret_cast<float*>(vB.buf[typeI])
                                [elemI] << " ";
                        }
                        else
                        {
                            vFile << reinterpret_cast<double*>(vB.buf[typeI])
                                [elemI] << " ";
                        }
                        elemI++;
                    }
                    vFile << std::endl;
                }
                for(label timeI = 0; timeI < vH.NbTimeStep; timeI++)
                {
                    for(label pointI = 0; pointI < nVertices_[typeI / 3];
                        pointI++)
                    {
                        for(label componentI = 0;
                            componentI < nComponents_[typeI % 3]; componentI++)
                        {
                            if(pF.size == 4)
                            {
                                vFile << reinterpret_cast<float*>
                                    (vB.buf[typeI])[elemI] << " ";
                            }
                            else
                            {
                                vFile << reinterpret_cast<double*>
                                    (vB.buf[typeI])[elemI] << " ";
                            }
                            elemI++;
                        }
                        vFile << std::endl;
                    }
                }
            }
        }

        // 2D texts
        for(label textI = 0, elemI = 0; textI < vH.NbT2; textI++)
        {
            for(label numI = 0; numI < 4; numI++)
            {
                if(pF.size == 4)
                {
                    vFile << reinterpret_cast<float*>(vB.T2D)[elemI] << " ";
                }
                else
                {
                    vFile << reinterpret_cast<double*>(vB.T2D)[elemI] << " ";
                }
                elemI++;
            }
            vFile << std::endl;
        }
        for(label charI = 0; charI < vH.t2l; charI++)
        {
            vFile << vB.T2C[charI];
        }
        vFile << std::endl;

        // 3D texts
        for(label textI = 0, elemI = 0; textI < vH.NbT3; textI++)
        {
            for(label numI = 0; numI < 5; numI++)
            {
                if(pF.size == 4)
                {
                    vFile << reinterpret_cast<float*>(vB.T3D)[elemI] << " ";
                }
                else
                {
                    vFile << reinterpret_cast<double*>(vB.T3D)[elemI] << " ";
                }
                elemI++;
            }
            vFile << std::endl;
        }
        for(label charI = 0; charI < vH.t3l; charI++)
        {
            vFile << vB.T3C[charI];
        }
    }

    vFile << std::endl;
    vFile << "$EndView" << std::endl;

    gInfo(verbosity_ >= 4) << "        ... done." << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
