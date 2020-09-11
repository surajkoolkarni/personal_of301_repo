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
    calcNusselt

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

        IOobject Theader
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check T exists
        if (Theader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading T" << endl;
            volScalarField T(Theader, mesh);

            Info<< " Calculating NusseltNumber" << endl;

            volScalarField NusseltNumber
            (
              IOobject
              (
                 "NusseltNumber",
                 runTime.timeName(),
                 mesh,                    
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
              ),  
            mesh,
            dimensionedScalar
            (
                 "NusseltNumber",
                 dimLength/dimLength,
                 0
            )
           );
        
           forAll(NusseltNumber.boundaryField(), patchi)
           {
            NusseltNumber.boundaryField()[patchi] = 
             length*mag((-T.boundaryField()[patchi].snGrad())/(T.boundaryField()[patchi]-T_amb+smallVal));
           }
    
           NusseltNumber.write();

            label patchID = mesh.boundaryMesh().findPatchID(patchName);

            const polyPatch& cPatch = mesh.boundaryMesh()[patchID];

            const surfaceScalarField& magSf = mesh.magSf();

            scalar Nu_sum = 0.0;
            scalar patchArea = 0.0;
    
             forAll(cPatch, facei)
             {
                Nu_sum += (magSf.boundaryField()[patchID][facei])*(NusseltNumber.boundaryField()[patchID][facei]);
                patchArea += (magSf.boundaryField()[patchID][facei]);
             }

             const double& NuAverage = Nu_sum/patchArea;

             Info << "    Average Nusselt Number is: " << NuAverage << endl; 
        
        } 
    
        else
        {
            Info<< "    No T" << endl;
        }
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
