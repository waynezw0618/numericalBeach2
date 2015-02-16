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

\*---------------------------------------------------------------------------*/

#include "numericalBeach.H"
#include "numericalBeachFvPatchField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

numericalBeach::numericalBeach(volVectorField& U)
:
    U_(U),
    
    IOdictionary
    (
     IOobject
     (
      "waveProperties",
      U.time().constant(),
      U.db(),
      IOobject::MUST_READ,
      IOobject::NO_WRITE
      )
     ),
    waveT_(lookup("waveT")),
    waveL_(lookup("waveL")),
    waveOmega_(2*mathematicalConstant::pi/waveT_),
    waveK_(2*mathematicalConstant::pi/waveL_),
    waveC_(waveL_/waveT_),
    waveTsoft_(lookupOrDefault<scalar>("waveTsoft",scalar(0.0))),
    massSource_(lookupOrDefault<Switch>("massSource",true)), //false
    momentumSource_(lookupOrDefault<Switch>("momentumSource",false)), //false
    secondOrder_(lookupOrDefault<Switch>("SecondOrder",false)),
    X_(lookupOrDefault<vector>("streamDirection",vector(1, 0, 0))),
    Z_(lookupOrDefault<vector>("streamDirection",vector(0, 1, 0))),
    Tdamp_(lookupOrDefault<tensor>("beachDirection",tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)-X_*X_)),
    h1_(lookupOrDefault<scalar>("h1",scalar(2.0))), // h1=0
    h2_(lookupOrDefault<scalar>("h2",scalar(0.0))), // h2> D+H/2 to include all cells in the vertical direction
    deltaL_(lookupOrDefault<scalar>("detL",scalar(0.016666))), // default value based on 60 cells/wave length
    Lx_(deltaL_* waveL_), //according to Gauss Theory, Source =2U/Lx,
    //zoneCenterXcoor_(lookupOrDefault<dimensionedScalar>("zoneCenterXcoor",dimensionedScalar("Xc",dimLength,0.0))),
   // AV_(lookupOrDefault<dimensionedScalar>("AV",dimensionedScalar("cctv",dimLength,0.0))),
    //AV_(zoneCenterXcoor_),
    zoneCenterZcoor_((h1_+h2_)/2*waveL_),
  //  curTsoft_(lookupOrDefault<scalar>("curTsoft",0.0)),
    sourceZone_(
                IOobject
                (
                 "sourceZone",
                 U.db().time().timeName(),
                 U.mesh(),
                 IOobject::NO_READ,
                 IOobject::NO_WRITE
                 ),
                U.mesh(),
                dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0), 0.0)
                )
{
  include "testgit.H"
}


numericalBeach::~numericalBeach()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> numericalBeach::damping() const
{
    tmp<volScalarField> tDamping
    (
        new volScalarField
        (
            IOobject
            (
                "damping",
                U_.db().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_.mesh(),
            dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0), 0.0)
        )
    );

    volScalarField& damping = tDamping();

    const volVectorField::GeometricBoundaryField& bvf = U_.boundaryField();
    forAll(bvf, patchi)
    {
        if (isType<numericalBeachFvPatchField>(bvf[patchi]))
        {
            const numericalBeachFvPatchField& beach =
                dynamic_cast<const numericalBeachFvPatchField&> (bvf[patchi]);

            damping = max(damping, beach.internalDamping());
        }
    }

    return tDamping;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
