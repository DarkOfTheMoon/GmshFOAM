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
    The class makes it possible to share text message lines between
    OpenFOAM applications and Gmsh by separately implimenting the
    stream writer function endl() for the message stream.

\*---------------------------------------------------------------------------*/

#include "gmshMessageStream.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * Static data * * * * * * * * * * * * * * //

bool gmshMessageStream::isInitialized_ = false;
gmshMessageStream gInfo(gmshMessageStream::gINFO);
gmshMessageStream gWarning(gmshMessageStream::gWARNING);
gmshMessageStream gSeriousError(gmshMessageStream::gSERIOUS);
gmshMessageStream gNull(gmshMessageStream::gNULL);

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

gmshMessageStream::gmshMessageStream(errorSeverity sev)
    : severity_(sev)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

gmshMessageStream& endl(gmshMessageStream& os)
{
    // we need to match the gmshMessageStream and the messageStream
    // severities here
    switch(os.severity_)
    {
        case gmshMessageStream::gINFO:
            Info << os.str_.str().c_str() << endl;
            break;
        case gmshMessageStream::gWARNING:
            Warning << os.str_.str().c_str() << endl;
            break;
        case gmshMessageStream::gSERIOUS:
            SeriousError << os.str_.str().c_str() << endl;
            break;
        case gmshMessageStream::gFATAL:
            FatalError << os.str_.str().c_str() << endl;
            break;
        case gmshMessageStream::gNULL:
            return os;
            break;
      }

    os.str_.str(std::string(""));

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
