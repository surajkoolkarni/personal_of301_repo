// - Added code here. -
   int iNormalDir(-1);
   word normalDirection;

   IOdictionary transportPropertiesDict
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
   transportPropertiesDict.lookup("normalDirection") >> normalDirection;

   if (normalDirection == "x")
       iNormalDir = 0;
   else if (normalDirection == "y")
       iNormalDir = 1;
   else if (normalDirection == "z")
       iNormalDir = 2;
   else
       FatalErrorIn(args.executable().c_str())
           << "Unknown Cartesian direction " << normalDirection << " attempting to be disabled!"
           << exit(FatalError);

   Vector<label>& solD = const_cast<Vector<label>&>(mesh.solutionD());
   Info << "Disabling solution in the " << normalDirection << " direction" << endl;
   solD[iNormalDir] = -1;
