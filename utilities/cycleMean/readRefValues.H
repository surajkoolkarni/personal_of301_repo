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

scalar nDirectories(readScalar(refValues.lookup("nDirectories")));
Info << "No of directories selected for averaging: "<< nDirectories << endl;
