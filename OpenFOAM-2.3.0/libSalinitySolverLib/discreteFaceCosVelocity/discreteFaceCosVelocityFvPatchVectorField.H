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

Class
    Foam::discreteFaceCosVelocityFvPatchVectorField

Description
    
    

SourceFiles
    discreteFaceCosVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef discreteFaceCosVelocityFvPatchVectorField_H
#define discreteFaceCosVelocityFvPatchVectorField_H

#include "mixedFvPatchFields.H"
#include "coordinateSystem.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class advectiveFvPatch Declaration
\*---------------------------------------------------------------------------*/

class discreteFaceCosVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
protected:

    // Private data

	//-
	label nIntervals_;

        //- 
        vector A_;
        
        //-
        scalar H_;
        
        //-
        vector Hdirection_;

        //- 
        scalar omega0_;
        
        //-
        scalar phi0_;
        
        //-
        scalar minZ_;

public:

    //- Runtime type information
    TypeName("discreteFaceCosVelocity");


    // Constructors

        //- Construct from patch and internal field
        discreteFaceCosVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        discreteFaceCosVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given discreteFaceCosVelocityFvPatchVectorField
        //  onto a new patch
        discreteFaceCosVelocityFvPatchVectorField
        (
            const discreteFaceCosVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        discreteFaceCosVelocityFvPatchVectorField
        (
            const discreteFaceCosVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField > clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new discreteFaceCosVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        discreteFaceCosVelocityFvPatchVectorField
        (
            const discreteFaceCosVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new discreteFaceCosVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access
        
        //-
        label whichInterval (scalar zfc);
        
        //-
        void setPatchVelocities(vectorField& pV);

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
        
        //-
        virtual void autoMap (const fvPatchFieldMapper& );
        
        //-
        virtual void rmap (const fvPatchVectorField& pf, const labelList& ll);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
