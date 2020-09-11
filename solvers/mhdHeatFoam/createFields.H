Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rho
(
    "rho",
    dimDensity,
    transportProperties
);

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties
);

dimensionedScalar mu
(
    "mu",
    dimensionSet(1, 1, -2, 0, 0, -2, 0),
    transportProperties
);

dimensionedScalar sigma
(
    "sigma",
    dimensionSet(-1, -3, 3, 0, 0, 2, 0),
    transportProperties
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

dimensionedScalar gammaT
(
    "gammaT",
    dimViscosity,
    transportProperties.lookup("gammaT")
);

#include "createPhi.H"

Info<< "Reading field pB\n" << endl;
volScalarField pB
(
    IOobject
    (
        "pB",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field B\n" << endl;
volVectorField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

// Info<<"Creating turbulence model\n"<<endl;
// autoPtr<incompressible::turbulenceModel> turbulence
// (
//     incompressible::turbulenceModel::New(U, phi, laminarTransport)
// );

#include "createPhiB.H"
// #include "readTransportProperties.H"
#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

dimensionedScalar DB = 1.0/(mu*sigma);
DB.name() = "DB";

dimensionedScalar DBU = 1.0/(2.0*mu*rho);
DBU.name() = "DBU";


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);

mesh.setFluxRequired(p.name());
mesh.setFluxRequired(pB.name());
