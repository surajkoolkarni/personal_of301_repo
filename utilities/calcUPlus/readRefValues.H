Info << "\nReading refValues" << endl;
IOdictionary refValues
(
    IOobject
        (
         "refValues",
         runTime.constant(),
         mesh,
         IOobject::MUST_READ,
         IOobject::NO_WRITE
        )
);



scalar nu(readScalar(refValues.lookup("nu")));
Info << "Kinematic viscosity is "<< nu << endl;
