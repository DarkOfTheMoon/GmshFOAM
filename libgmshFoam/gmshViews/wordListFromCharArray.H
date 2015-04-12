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
    wordListFromCharArray

Description
    creates wordList from an array of (char*)'s.

SourceFiles
    (none)

\*---------------------------------------------------------------------------*/

#ifndef wordListFromCharArray_H
#define wordListFromCharArray_H

#include "wordList.H"

namespace Foam
{

class wordListFromCharArray
: public wordList
{
public:
    // constructor
    wordListFromCharArray(const label size, const char *cA[]) : wordList()
    {
        setSize(size);
        forAll(*this, wordI)
        {
            operator[](wordI) = cA[wordI];
        }
    }            
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
