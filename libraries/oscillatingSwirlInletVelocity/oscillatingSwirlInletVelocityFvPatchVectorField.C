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

#include "oscillatingSwirlInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oscillatingSwirlInletVelocityFvPatchVectorField::
oscillatingSwirlInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    rpm_(),
    freq_(),
    Umax_(),
    cavityRadius_(vector::zero),
    y_(vector::zero)
{}


Foam::oscillatingSwirlInletVelocityFvPatchVectorField::
oscillatingSwirlInletVelocityFvPatchVectorField
(
    const oscillatingSwirlInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    rpm_(ptf.rpm_().clone().ptr()),
    freq_(ptf.freq_().clone().ptr()),
    Umax_(ptf.Umax_().clone().ptr()),
    cavityRadius_(ptf.cavityRadius_().clone().ptr()),
    y_(ptf.y_().clone().ptr)
{}


Foam::oscillatingSwirlInletVelocityFvPatchVectorField::
oscillatingSwirlInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    rpm_(DataEntry<scalar>::New("rpm", dict)),
    freq_(DataEntry<scalar>::New("freq",dict)),
    Umax_(DataEntry<scalar>::New("Umax",dict)),
    cavityRadius_(DataEntry<vector>::New("cavityRadius",dict)),
    y_(DataEntry<vector>::New("y",dict))
{}


Foam::oscillatingSwirlInletVelocityFvPatchVectorField::
oscillatingSwirlInletVelocityFvPatchVectorField
(
    const oscillatingSwirlInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    rpm_(ptf.rpm_().clone().ptr()),
    freq_(ptf.freq_().clone().ptr()),
    Umax_(ptf.Umax_().clone().ptr()),
    cavityRadius_(ptf.cavityRadius_().clone().ptr()),
    y_(ptf.y_().clone().ptr)
{}


Foam::oscillatingSwirlInletVelocityFvPatchVectorField::
oscillatingSwirlInletVelocityFvPatchVectorField
(
    const oscillatingSwirlInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    rpm_(ptf.rpm_().clone().ptr()),
    freq_(ptf.freq_().clone().ptr()),
    Umax_(ptf.Umax_().clone().ptr()),
    cavityRadius_(ptf.cavityRadius_().clone().ptr()),
    y_(ptf.y_().clone().ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::oscillatingSwirlInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    boundBox bb(patch().patch().localPoints(), true);

    const scalar t = this->db().time().timeOutputValue();
    const scalar rpm = rpm_->value(t);
    const scalar freq = freq_->value(t);

    const scalar totArea = gSum(patch().magSf());

    const vector avgCenter = gSum(patch().Cf()*patch().magSf())/totArea;
    const vector avgNormal = gSum(patch().Sf())/totArea;

    const scalar Umax = Umax_->value(t);
    const scalar sineval = sin(2*(constant::mathematical::pi)*freq*t);

    const vectorField& c = patch().Cf();

    scalarField coord = (c & y_)/((bb.max()) & y_);


    // Update angular velocity - convert [rpm] to [rad/s]
    tmp<vectorField> tangentialVelocity
        (
            (rpm*constant::mathematical::pi/30.0)
          * (patch().Cf() - avgCenter) ^ avgNormal
        );

    tmp<vectorField> n = patch().nf();

    if (sineval>=0)
    {
        // swirl + pulsations during expulsion phase
        operator==(tangentialVelocity + n*Umax*sineval*(1-sqr(coord)));
    }
    else
    {
        // no swirl during suction phase
         operator==((n*Umax*sineval)*(1-sqr(coord)));
    }

    // else
    // {
    //     FatalErrorIn
    //     (
    //         "oscillatingSwirlInletVelocityFvPatchVectorField::updateCoeffs()"
    //     )
    //         << "    on patch " << this->patch().name()
    //         << " of field " << this->dimensionedInternalField().name()
    //         << " in file " << this->dimensionedInternalField().objectPath()
    //         << nl << exit(FatalError);
    // }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::oscillatingSwirlInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);
    rpm_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       oscillatingSwirlInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
