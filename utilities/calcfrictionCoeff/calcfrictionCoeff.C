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
    calcCf

Description
    Calculates and writes skin friction coefficient without including turbulent
    shear stress.

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

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject Uxheader
        (
            "Ux",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject Uyheader
        (
            "Uy",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check U exists
        if (Uxheader.headerOk())
        {
          if(Uyheader.headerOk())
          {
            mesh.readUpdate();

            Info<< "    Reading U" << endl;
            volScalarField Ux(Uxheader, mesh);
            volScalarField Uy(Uyheader, mesh);

            Info<< " Calculating frictionCoeff" << endl;

            volScalarField frictionCoeff
            (
              IOobject
              (
                 "frictionCoeff",
                 runTime.timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
              ),
            mesh,
            dimensionedScalar
            (
                 "frictionCoeff",
                 dimLength/dimLength,
                 0
            )
           );

           volScalarField gradUx
           (
             IOobject
             (
                "gradUx",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
             ),
           mesh,
           dimensionedScalar
           (
                "gradUx",
                dimVelocity/dimLength,
                0
           )
           );

           volScalarField gradUy
           (
             IOobject
             (
                "gradUy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
             ),
           mesh,
           dimensionedScalar
           (
                "gradUy",
                dimVelocity/dimLength,
                0
           )
           );

           forAll(gradUx.boundaryField(), patchi)
           {
             gradUx.boundaryField()[patchi] = Ux.boundaryField()[patchi].snGrad();
             frictionCoeff.boundaryField()[patchi] =
             nu_kine*(sqrt(pow((-Ux.boundaryField()[patchi].snGrad()),2)+pow((-Uy.boundaryField()[patchi].snGrad()),2)))/dynamicPressure;
           }
           forAll(gradUx.boundaryField(), patchi)
           {
             gradUy.boundaryField()[patchi] = Uy.boundaryField()[patchi].snGrad();
           }
           frictionCoeff.write();
           gradUx.write();
           gradUy.write();

            label patchID = mesh.boundaryMesh().findPatchID(patchName);

            const polyPatch& cPatch = mesh.boundaryMesh()[patchID];

            const surfaceScalarField& magSf = mesh.magSf();

            scalar Cf_sum = 0.0;
            scalar patchArea = 0.0;

             forAll(cPatch, facei)
             {
                Cf_sum += (magSf.boundaryField()[patchID][facei])*(frictionCoeff.boundaryField()[patchID][facei]);
                patchArea += (magSf.boundaryField()[patchID][facei]);
             }

             const double& CfAverage = Cf_sum/patchArea;

             Info << "    Average frictionCoeff is: " << CfAverage << endl;
           }
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
