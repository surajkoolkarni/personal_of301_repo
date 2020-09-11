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

scalar T_amb(readScalar(refValues.lookup("T_amb")));
Info << "Ambient temperature is:"<< T_amb << endl;

scalar length(readScalar(refValues.lookup("length")));
Info << "Length scale is set to:"<< length << endl;

word patchName(word(refValues.lookup("patchName"))); 