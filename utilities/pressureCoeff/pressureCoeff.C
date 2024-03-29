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
    pressureCoeff

Description
    Calculates and writes the volVectorField of squared pressure.

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

        IOobject pheader
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (pheader.headerOk())
        {
            if (Uheader.headerOk())
            {
              mesh.readUpdate();

              Info<< "    Reading p" << endl;
              volScalarField p(pheader, mesh);

              Info<< "    Reading U" << endl;
              volVectorField U(Uheader,mesh);

              Info<< " Calculating pressureCoeff   " << endl;

              volScalarField pressureCoeff
              (
                IOobject
                (
                   "pressureCoeff",
                   runTime.timeName(),
                   mesh,
                   IOobject::NO_READ,
                   IOobject::AUTO_WRITE
                ),
              mesh,
              dimensionedScalar
              (
                   "pressureCoeff",
                   dimensionSet(0,0,0,0,0,0,0),
                   0
              )
              );

             forAll(pressureCoeff.boundaryField(),patchI)
             {
               pressureCoeff.boundaryField()[patchI] =
                (p.boundaryField()[patchI]- pRef)/(smallVal+Uimp*Uimp);
             }

              pressureCoeff.write();
            }

            else
            {
              Info<<"     No U"<<endl;
            }
        }

        else
        {
            Info<< "    No p" << endl;
        }
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
