/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C)
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
    Monitors a OpenFOAM case directory and informs Gmsh everytime new
    time directory is created.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#ifdef cygwin
#include "Time.hh"
#else
#include "Time.H"
#endif

#include "GmshClient.h"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef enum { sendOptions, monitor } action;

bool readDictionaryOptions(string& viewOptionStr, const string& optionName,
const Time& runTime)
{
    string optionStr;

    objectRegistry objRegistry(runTime);
    IOdictionary monitorDict
    (
        IOobject("gmshFoamMonitorDict", runTime.system(), objRegistry,
        IOobject::READ_IF_PRESENT, IOobject::NO_WRITE)
    );
    if(monitorDict.headerOk() && monitorDict.found(optionName))
    {
        viewOptionStr = monitorDict.lookup(optionName);
        return true;
    }
    return false;
}

void doMonitoring(const fileName& rootPath, const fileName& caseName,
GmshClient& client, const label interval)
{
    Time runTime(Time::controlDictName, rootPath, caseName);

    fileName casePath = rootPath/caseName;

    scalar stopAt;

    string stopAtStr;
    if(readDictionaryOptions(stopAtStr, "stopAt", runTime))
    {
        if(stopAtStr == "endTime")
        {
            stopAt = runTime.endTime().value();
        }
        else
        {
            stopAt = atof(stopAtStr.c_str());
        }
    }
    else
    {
        stopAt = runTime.endTime().value();
    }

    OStringStream infoStr;
    infoStr << "Starting monitoring the case directory " << casePath.c_str()
        << " up to stopAt = " << stopAt << " with monitoring interval "
        << interval << " seconds." << endl;
    client.Info(infoStr.str().c_str());

    OStringStream initializeStr;
    initializeStr << "View.ShowTime=2;" << endl
        << "PostProcessing.FoamStartTime=\"latestTime\";" << endl;

    string initialOptionStr;
    if(readDictionaryOptions(initialOptionStr, "initialViewOptions", runTime))
    {
        initializeStr << initialOptionStr.c_str() << endl;
    }

    client.ParseString(initializeStr.str().c_str());

    label timeI = 0;
    while(1)
    {
        instantList timeList(runTime.times());

        if(timeI < timeList.size())
        {
            const instant& latestInstant = timeList[timeList.size() - 1];

            // notify user of a new field data
            time_t t = time(NULL);
            OStringStream newDirectoryStr;
            newDirectoryStr << ctime(&t) << ": "
                << "A new directory t = "
                << latestInstant.name().c_str() << " found." << endl;
            client.Info(newDirectoryStr.str().c_str());

            // instruct Gmsh to read the newest time directory
            OStringStream mergeStr;
            mergeStr << "For viewI In {0:PostProcessing.NbViews-1}" << endl
                << "Delete View[0];" << endl
                << "EndFor" << endl;
            mergeStr << "Merge \"" << casePath.c_str()
                << "/system/controlDict\";" << endl;

            // read the view option string everytime the views are
            // renewed so that the new options can be applied
            // on-the-fly
            string viewOptionStr;
            if(readDictionaryOptions(viewOptionStr, "viewOptions", runTime))
            {
                mergeStr << viewOptionStr.c_str() << endl;
            }

            mergeStr << "Draw;" << endl;

            client.ParseString(mergeStr.str().c_str());

            if(latestInstant.value() >= stopAt - SMALL)
            {
                // finish and clean-up
                OStringStream reachedEndStr;
                reachedEndStr << "Reached the stopAt t = " << stopAt
                    << ". Finishing monitoring." << endl;
                client.Info(reachedEndStr.str().c_str());

                client.Stop();
                client.Disconnect();
                return;
            }

            timeI = timeList.size();
        }

        // sleep for the specified interval
        Foam::sleep(interval);
    }
}

// Main program:
int main(int argc, char *argv[])
{
    argList::validOptions.insert("monitor", "interval");
    argList::validOptions.insert("socket", "socket");

#include "setRootCase.H"

    action whatToDo = sendOptions;
    string socket, intervalStr;

    if(args.options().found("monitor"))
    {
        intervalStr = args.options()["monitor"];
        whatToDo = monitor;
    }
    if(args.options().found("socket"))
    {
        socket = args.options()["socket"];
    }
    else
    {
        socket = fileName(getenv("HOME"))/".gmshsock";
    }

    GmshClient client;
    if(client.Connect(socket.c_str()) < 0)
    {
        FatalErrorIn(args.executable())
            << "Unable to connect to Gmsh socket " << socket.c_str() << "."
                << exit(FatalError);
    }
    client.Start();

    if(whatToDo == sendOptions)
    {
        // selection of monitoring interval in seconds
        client.Option(1, "5");
        client.Option(1, "15");
        client.Option(1, "60");
    }
    else if(whatToDo == monitor)
    {
        // do monitoring!
        const label interval = atoi(intervalStr.c_str());
        doMonitoring(args.rootPath(), args.caseName(), client, interval);
    }

    return 0;
}

// ************************************************************************* //
