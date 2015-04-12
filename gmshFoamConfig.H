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

Description
    Simple compatibility mechanism across different versions of OpenFOAM.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef gmshFoamConfig_H
#define gmshFoamConfig_H

// OF >= 1.4 have Time(root, case) constructor which doesn't require
// existence of controlDict and symmTensor
#if WM_PROJECT_VERSION_MAJOR > 1 || (WM_PROJECT_VERSION_MAJOR == 1 \
    && WM_PROJECT_VERSION_MINOR >= 4)
#define WITH_NOREAD_TIME_CTOR 1
#define WITH_SYMMTENSOR 1
#else
#define WITH_NOREAD_TIME_CTOR 0
#define WITH_SYMMTENSOR 0
#endif

// OF >= 1.4.1 have fvsPatchField
#if WM_PROJECT_VERSION_MAJOR == 1 && (WM_PROJECT_VERSION_MINOR < 4 \
    || WM_PROJECT_VERSION_MINOR == 4 && WM_PROJECT_VERSION_PATCH == 0)
#define WITH_FVSPATCHFIELD 0
#else
#define WITH_FVSPATCHFIELD 1
#endif

#endif
