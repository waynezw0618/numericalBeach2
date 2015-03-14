/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    numerical damping in the internalfield for selected boundaries

Class
    Foam::numericalBeach

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved.

SourceFiles
    numericalBeach.C

\*---------------------------------------------------------------------------*/

#ifndef numericalBeach_H
#define numericalBeach_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class numericalBeach Declaration
\*---------------------------------------------------------------------------*/

class numericalBeach
:
    public IOdictionary
{

    // Private data

    //- Velocity field to be damped at boundaries
    volVectorField& U_;
    //- for test only
    //- wave parameters
   /* dimensionedScalar waveT_;
    dimensionedScalar waveL_;
    dimensionedScalar waveOmega_;
    dimensionedScalar waveK_;
    dimensionedScalar waveC_;
    scalar waveTsoft_;
    
    Switch massSource_;
    Switch momentumSource_;
    Switch secondOrder_;
    
    // Read geometry parameters
    vector X_;
    vector Z_;
    tensor Tdamp_;
    
    scalar h1_; // h1=0
    scalar h2_; // h2> D+H/2 to include all cells in the vertical direction
    scalar deltaL_;
    dimensionedScalar Lx_;*/
    dimensionedScalar zoneCenterXcoor_;
  //  dimensionedScalar zoneCenterZcoor_;
    
    // for Beach current
   // dimensionedScalar AV_;
  //  scalar curTsoft_;

//    volScalarField sourceZone_;





    // Private Member Functions

        //- Disallow default bitwise copy construct
        numericalBeach(const numericalBeach&);

        //- Disallow default bitwise assignment
        void operator=(const numericalBeach&);


public:

    // Constructors

        //- Construct from components
        numericalBeach(volVectorField& U);


    // Destructor

        ~numericalBeach();


    // Member Functions

        // Access

        //- Return the source term for the momentum equation
        tmp<volScalarField> damping() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
