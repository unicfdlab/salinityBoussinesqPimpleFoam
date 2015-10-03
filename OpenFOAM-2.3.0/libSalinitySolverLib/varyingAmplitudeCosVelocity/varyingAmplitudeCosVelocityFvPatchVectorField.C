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

#include "varyingAmplitudeCosVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::varyingAmplitudeCosVelocityFvPatchVectorField::varyingAmplitudeCosVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    At_(),
    H_(1.0),
    Hdirection_(0.0, 0.0, 0.0),
    omega0_(0.0),
    phi0_(0.0),
    minZ_(0.0)
{
    this->refValue() = pTraits<vector>::zero;
    this->refGrad() = pTraits<vector>::zero;
    this->valueFraction() = 0.0;
    this->setPatchVelocities(this->refValue());
    
    Info << "Executing varyingAmplitudeCosVelocityFvPatchVectorField(fvPatch, DimensionedField) " << endl;
}


Foam::varyingAmplitudeCosVelocityFvPatchVectorField::varyingAmplitudeCosVelocityFvPatchVectorField
(
    const varyingAmplitudeCosVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<vector>(ptf, p, iF, mapper),
    At_(ptf.At_),
    H_(ptf.H_),
    Hdirection_(ptf.Hdirection_),
    omega0_(ptf.omega0_),
    phi0_(ptf.phi0_),
    minZ_(ptf.minZ_)
{
    this->setPatchVelocities(this->refValue());
    mixedFvPatchVectorField::updateCoeffs();

    Info << "Executing varyingAmplitudeCosVelocityFvPatchVectorField(varyingAmplitudeCosVelocity, fvPatch, DimensionedField, fvPatchFieldMapper) " << endl;
}


Foam::varyingAmplitudeCosVelocityFvPatchVectorField::varyingAmplitudeCosVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<vector>(p, iF),
    At_(dict)
{
    dict.lookup("H") >> H_;
    dict.lookup("Hdirection") >> Hdirection_;
    dict.lookup("omega0") >> omega0_;
    dict.lookup("phi0") >> phi0_;
    dict.lookup("minZ") >> minZ_;

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            Field<vector>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<vector>::operator=(this->patchInternalField());
    }

    this->refValue() = pTraits<vector>::zero;
    this->refGrad() = pTraits<vector>::zero;
    this->valueFraction() = 1.0;
    this->setPatchVelocities(this->refValue());
    mixedFvPatchVectorField::updateCoeffs();

    Info << "Executing varyingAmplitudeCosVelocityFvPatchVectorField(fvPatch, DimensionedField, dictionary) " << endl;
}


Foam::varyingAmplitudeCosVelocityFvPatchVectorField::varyingAmplitudeCosVelocityFvPatchVectorField
(
    const varyingAmplitudeCosVelocityFvPatchVectorField& ptpsf
)
:
    mixedFvPatchVectorField(ptpsf),
    At_(ptpsf.At_),
    H_(ptpsf.H_),
    Hdirection_(ptpsf.Hdirection_),
    omega0_(ptpsf.omega0_),
    phi0_(ptpsf.phi0_),
    minZ_(ptpsf.minZ_)
{
    this->setPatchVelocities(this->refValue());
    Info << "Executing varyingAmplitudeCosVelocityFvPatchVectorField(const varyingAmplitudeCosVelocity&) " << endl;
}


Foam::varyingAmplitudeCosVelocityFvPatchVectorField::varyingAmplitudeCosVelocityFvPatchVectorField
(
    const varyingAmplitudeCosVelocityFvPatchVectorField& ptpsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(ptpsf, iF),
    At_(ptpsf.At_),
    H_(ptpsf.H_),
    Hdirection_(ptpsf.Hdirection_),
    omega0_(ptpsf.omega0_),
    phi0_(ptpsf.phi0_),
    minZ_(ptpsf.minZ_)
{
    this->setPatchVelocities(this->refValue());

    Info << "Executing varyingAmplitudeCosVelocityFvPatchVectorField(const varyingAmplitudeCosVelocity, const DimensionedField) " << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::varyingAmplitudeCosVelocityFvPatchVectorField::autoMap(const fvPatchFieldMapper& m)
{
    mixedFvPatchField<vector>::autoMap(m);

    Info << "Executing autoMap " << endl;
}

void Foam::varyingAmplitudeCosVelocityFvPatchVectorField::rmap(const fvPatchVectorField& pf, const labelList& ll)
{
    mixedFvPatchField<vector>::rmap(pf, ll);
    Info << "Executing rmap " << endl;
}

void Foam::varyingAmplitudeCosVelocityFvPatchVectorField::setPatchVelocities (vectorField& pV)
{
    scalarField z = (this->patch().Cf() & Hdirection_) - minZ_;
    
    scalar t = this->patch().boundaryMesh().mesh().time().value();
    vector cA = At_(t);
    
    forAll(pV, iFace)
    {
	this->refGrad()[iFace] = vector (0.0, 0.0, 0.0);
	this->valueFraction()[iFace] = 1.0;
	pV[iFace] = 
	    cA * cos(Foam::constant::mathematical::pi * z[iFace] / H_) * (-omega0_) * sin(omega0_ * t + phi0_);
	if (z[iFace] < 0.0 || z[iFace] > H_)
	{
	    pV[iFace] = cA*0.0;
	}
    }
}

void Foam::varyingAmplitudeCosVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
    this->setPatchVelocities(this->refValue());
    
    Info << "max/min velocity: " << gMax(this->refValue()) << "/" << gMin(this->refValue()) << endl;

    mixedFvPatchField<vector>::updateCoeffs();
}


void Foam::varyingAmplitudeCosVelocityFvPatchVectorField::write(Ostream& os) const
{

    mixedFvPatchVectorField::write(os);

    At_.write(os);
        
    os.writeKeyword("H") << H_ << token::END_STATEMENT << nl;
    
    os.writeKeyword("minZ") << minZ_ << token::END_STATEMENT << nl;
    
    os.writeKeyword("Hdirection") << Hdirection_ << token::END_STATEMENT << nl;
    
    os.writeKeyword("omega0") << omega0_ << token::END_STATEMENT << nl;
    
    os.writeKeyword("phi0") << phi0_ << token::END_STATEMENT << nl;
    
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        varyingAmplitudeCosVelocityFvPatchVectorField
    );
}



// ************************************************************************* //
