/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    U2

Description
    Calculates and writes the volVectorField of squared velocity.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject Uyheader
        (
            "Uy",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check Uy exists
        if (Uyheader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading Uy" << endl;
            volScalarField Uy(Uyheader, mesh);

            Info<< " Calculating Squared velocity U2y   " << endl;

            volScalarField U2y
            (
              IOobject
              (
                 "U2y",
                 runTime.timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
              ),
            mesh,
            dimensionedScalar
            (
                 "U2y",
                 dimVelocity*dimVelocity,
                 0
            )
            );

         forAll(U2y.internalField(),cellI)
         {
             U2y.internalField()[cellI] =
              (Uy.internalField()[cellI])*(Uy.internalField()[cellI]);
         }

          U2y.write();

        }

        else
        {
            Info<< "    No U" << endl;
        }
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
