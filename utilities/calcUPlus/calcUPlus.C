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
    calcDeltaS

Description
    Calculates and writes the gradient of T and nusselt number at the
    boundaries.

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
    #include "readRefValues.H"

    const double smallVal = 1e-10;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check T exists
        if (Uheader.headerOk())
        {

          Info<< "    Reading U" << endl;
          volVectorField U(Uheader, mesh);

          volScalarField uTau
          (
            IOobject
            (
              "uTau",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
            ),
           mesh,
           dimensionedScalar
           (
             "uTau",
             dimVelocity,
             0
           )
          );

            volScalarField DeltaS
            (
              IOobject
              (
                 "DeltaS",
                 runTime.timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
              ),
            mesh,
            dimensionedScalar
            (
                 "DeltaS",
                 dimLength,
                 0
            )
           );

           Info<< " Calculating Friction velocity" << endl;

           forAll(uTau.boundaryField(),patchi)
           {
             uTau.boundaryField()[patchi] =
             sqrt((nu*(mag(U.boundaryField()[patchi].snGrad()))));
           }

           uTau.write();

           Info<< "  Calculating DeltaS" << endl;

           forAll(DeltaS.boundaryField(), patchi)
           {
            DeltaS.boundaryField()[patchi] =
             sqrt(nu/(smallVal+(mag((U.boundaryField()[patchi].snGrad())))));
           }

           DeltaS.write();

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
