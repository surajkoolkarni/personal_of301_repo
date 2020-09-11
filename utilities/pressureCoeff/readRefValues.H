Info << " \nReading refValues" << endl;
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

scalar pRef(readScalar(refValues.lookup("pRef")));
Info << "Ambient pressure is:"<< pRef << endl;

scalar Uimp(readScalar(refValues.lookup("Uimp")));
Info << "Impingement velocity is:"<< Uimp << endl;

word patchName(word(refValues.lookup("patchName")));
