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

scalar dynamicPressure(readScalar(refValues.lookup("dynamicPressure")));

scalar length(readScalar(refValues.lookup("length")));
Info << "Length scale is set to:"<< length << endl;

scalar nu_kine(readScalar(refValues.lookup("nu_kine")));
Info << "Kinematic viscosity is:" << nu_kine << endl;

word patchName(word(refValues.lookup("patchName")));
