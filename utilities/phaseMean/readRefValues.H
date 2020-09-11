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

scalar nSpan(readScalar(refValues.lookup("nSpan")));
Info << "No of directories selected for phase averaging: "<< nSpan << endl;

scalar nTotal(readScalar(refValues.lookup("nTotal")));
Info << "No of total time steps are: " << nTotal <<endl;