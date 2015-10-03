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

#include "buoyantEqFluxPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buoyantEqFluxPressureFvPatchScalarField::
buoyantEqFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    rhoName_("rho"),
    fluxEqName_("fluxEq")
{}


Foam::buoyantEqFluxPressureFvPatchScalarField::
buoyantEqFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    fluxEqName_(dict.lookupOrDefault<word>("fluxEq", "fluxEq"))
{
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
}


Foam::buoyantEqFluxPressureFvPatchScalarField::
buoyantEqFluxPressureFvPatchScalarField
(
    const buoyantEqFluxPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    rhoName_(ptf.rhoName_),
    fluxEqName_(ptf.fluxEqName_)
{}


Foam::buoyantEqFluxPressureFvPatchScalarField::
buoyantEqFluxPressureFvPatchScalarField
(
    const buoyantEqFluxPressureFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    rhoName_(ptf.rhoName_),
    fluxEqName_(ptf.fluxEqName_)
{}


Foam::buoyantEqFluxPressureFvPatchScalarField::
buoyantEqFluxPressureFvPatchScalarField
(
    const buoyantEqFluxPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    rhoName_(ptf.rhoName_),
    fluxEqName_(ptf.fluxEqName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::buoyantEqFluxPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const uniformDimensionedVectorField& g =
        db().lookupObject<uniformDimensionedVectorField>("g");

    const fvPatchField<scalar>& rho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    // If the variable name is "p_rgh", "ph_rgh" or "pd"
    // assume it is p? - rho*g.h and set the gradient appropriately.
    // Otherwise assume the variable is the static pressure.
    if
    (
        dimensionedInternalField().name() == "p_rgh"
     || dimensionedInternalField().name() == "ph_rgh"
     || dimensionedInternalField().name() == "pd"
    )
    {
        gradient() = -rho.snGrad()*((g.value() & patch().Cf()));
        //if (patch().boundaryMesh().mesh().thisDb().foundObject<surfaceScalarField>(fluxEqName_))
        {
	    const fvsPatchField<scalar> fluxEq = patch().lookupPatchField<surfaceScalarField, scalar> (fluxEqName_);
	    gradient() -= (fluxEq / patch().magSf());
        }
    }
    else
    {
        gradient() = rho*((g.value() & patch().nf()));
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::buoyantEqFluxPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "fluxEq", "fluxEq", fluxEqName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        buoyantEqFluxPressureFvPatchScalarField
    );
}


// ************************************************************************* //
