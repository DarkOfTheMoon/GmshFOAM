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

void gmshViewBase::getViewSize(gmshViewHeader& vH, label vSize[])
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
