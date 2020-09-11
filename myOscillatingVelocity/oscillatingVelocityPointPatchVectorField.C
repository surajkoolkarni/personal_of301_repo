/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "oscillatingVelocityPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oscillatingVelocityPointPatchVectorField::
oscillatingVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    amplitude_(0.0),
    omega_(0.0),
    rc_(0.0),
    p0_(p.localPoints())
{}


oscillatingVelocityPointPatchVectorField::
oscillatingVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    omega_(readScalar(dict.lookup("omega"))),
    rc_(readScalar(dict.lookup("rc")))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dict , p.size());
    }
    else
    {
        p0_ = p.localPoints();
    }
}


oscillatingVelocityPointPatchVectorField::
oscillatingVelocityPointPatchVectorField
(
    const oscillatingVelocityPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_),
    omega_(ptf.omega_),
    rc_(ptf.rc_),
    p0_(ptf.p0_, mapper)
{}


oscillatingVelocityPointPatchVectorField::
oscillatingVelocityPointPatchVectorField
(
    const oscillatingVelocityPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    amplitude_(ptf.amplitude_),
    omega_(ptf.omega_),
    rc_(ptf.rc_),
    p0_(ptf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void oscillatingVelocityPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    p0_.autoMap(m);
}


void oscillatingVelocityPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const oscillatingVelocityPointPatchVectorField& oVptf =
        refCast<const oscillatingVelocityPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(oVptf, addr);

    p0_.rmap(oVptf.p0_, addr);
}


void oscillatingVelocityPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->dimensionedInternalField().mesh()();
    const Time& t = mesh.time();
    const pointPatch& p = this->patch();
    vectorField sd = p0_;
    vector p0rot;
    scalar pointMove = 0.0;
    
        
    forAll(p0_,iter)
    {
        p0rot = p0_[iter]; // p relative to new origin rotated
        // Plane from x and y values calculated and inserted into z values
        
        //pointMove  = 0.003 - 5.926*p0rot[0]*p0rot[0]-5.926*p0rot[2]*p0rot[2];
        
        pointMove  = amplitude_*(1-(p0rot[0]*p0rot[0])/(rc_*rc_)-(p0rot[2]*p0rot[2])/(rc_*rc_))*sin(omega_*t.value());
        
        
        p0rot = vector(0,pointMove,0);//*sin(omega_*t.value());//*sin(2*(22/7)*1.25*t.value());
        

       // p0rot = vector(0, 0, p0rot[0]);//+X1_*p0rot[0]+Y2_*p0rot[1]*p0rot[1]+Y1_*p0rot[1]+Cconst_)*sin(25.12*defTime_*t.value());
        sd[iter] = p0rot ; // Plane rotated back to original position
        
        /*
        if(iter ==30)
        {
        Info<< "t.value() = " << t.value() << " sin = " << sin(7.85398*t.value()) << nl << endl;
        Info<< "iter = " << iter << " pointMove = " << pointMove << nl << endl;
        Info<< "sin*pointMove = " << pointMove << nl << endl;
        Info<< "sd[iter]  = " << sd[iter]  << nl << endl;
	    }*/
	    
    };
    
    Field<vector>::operator=
    (
    //sd * sin(2*(22/7)*defTime_*t.value())
   // p0_ + sd 
    
        (p0_ + sd - p.localPoints())
       /t.deltaTValue()
       
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void oscillatingVelocityPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeKeyword("amplitude")
        << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("omega")
        << omega_ << token::END_STATEMENT << nl;
    os.writeKeyword("rc")
        << rc_ << token::END_STATEMENT << nl;
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    oscillatingVelocityPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
