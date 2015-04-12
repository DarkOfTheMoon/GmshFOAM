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
    gmshMessageStream

Description
    The class makes it possible to share text message lines between
    OpenFOAM applications and Gmsh by separately implimenting the
    stream writer function endl() for the message stream.

SourceFiles
    gmshMessageStream.C (for OpenFOAM applications)
    gmsh/Geo/GModelIO_FOAM.cpp (for Gmsh)

\*---------------------------------------------------------------------------*/

#ifndef gmshMessageStream_H
#define gmshMessageStream_H

#include <sstream>
#include "OStringStream.H"
#include "messageStream.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class gmshMessageStream;

extern gmshMessageStream gNull;
extern gmshMessageStream gInfo;
extern gmshMessageStream gWarning;
extern gmshMessageStream gSeriousError;

gmshMessageStream& endl(gmshMessageStream& os);

// the text message stream writer class
class gmshMessageStream
{
public:
    enum errorSeverity
    {
        gINFO,       // Debugging information in event of error
        gWARNING,    // Warning of possible problem
        gSERIOUS,    // A serious problem (data corruption?)
        gFATAL,      // Oh bugger!
        gNULL        // Dummy stream
    };

private:
    // Private data

    // using std::stringstream because Foam::OStringStream doesn't have a
    // clearance method.
    std::stringstream str_;
    errorSeverity severity_;

    static bool isInitialized_;

public:
    gmshMessageStream(errorSeverity sev);

    // for invoking the endl() friend function
    gmshMessageStream& operator<<(gmshMessageStream& (*f)(gmshMessageStream &))
        { return f(*this); }
    // a super-easy way to make use of operator<<()'s defined for many Foam
    // classes
    template <class T> gmshMessageStream& operator<<(const T t)
        { if(severity_ != gNULL) { OStringStream os; os << t;
        str_ << os.str().c_str(); } return *this; }

    // Evaluate the given logical expression and returns the gNull
    // stream if false
    gmshMessageStream& operator()(const bool b)
        { return b ? *this : gNull; }

    // the stream writer; have different implimentations for gmsh2ToFoam and
    // Gmsh
    friend gmshMessageStream& endl(gmshMessageStream& os);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
