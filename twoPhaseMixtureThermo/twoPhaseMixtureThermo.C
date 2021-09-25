/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "twoPhaseMixtureThermo.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
#include "collatedFileOperation.H"

#include "liquidThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseMixtureThermo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseMixtureThermo::twoPhaseMixtureThermo
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    psiThermo(U.mesh(), word::null),
    twoPhaseMixture(U.mesh(), *this),
    interfaceProperties(alpha1(), U, *this),
    thermo1_(nullptr),
    thermo2_(nullptr),
    mesh_(U.mesh()),
    pressureX_(10000, 0.0),
    temperatureY_(10000, 0.0)
{
    {
        volScalarField T1
        (
            IOobject
            (
                IOobject::groupName("T", phase1Name()),
                U.mesh().time().timeName(),
                U.mesh()
            ),
            T_,
            calculatedFvPatchScalarField::typeName
        );
        T1.write();
    }

    {
        volScalarField T2
        (
            IOobject
            (
                IOobject::groupName("T", phase2Name()),
                U.mesh().time().timeName(),
                U.mesh()
            ),
            T_,
            calculatedFvPatchScalarField::typeName
        );
        T2.write();
    }

    // Note: we're writing files to be read in immediately afterwards.
    //       Avoid any thread-writing problems.
    fileHandler().flush();

    thermo1_ = rhoThermo::New(U.mesh(), phase1Name());

    thermo2_ = rhoThermo::New(U.mesh(), phase2Name());

    // thermo1_->validate(phase1Name(), "e");
    // thermo2_->validate(phase2Name(), "e");

    correct();

    constructVaporPressureTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseMixtureThermo::~twoPhaseMixtureThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::twoPhaseMixtureThermo::constructVaporPressureTable()
{
    // Temperature range: [Tmin, Tmax)
    const scalar Tmin = 100.0;
    const scalar Tmax = 500.0;
    const label length = temperatureY_.size();
    forAll (temperatureY_, i)
    {
        temperatureY_[i] = Tmin + (Tmax - Tmin)/length*i;
        pressureX_[i] = pv(0.0, temperatureY_[i]);
    }
}

Foam::scalar Foam::twoPhaseMixtureThermo::Tsat(const scalar& tp) const
{
    const label length = temperatureY_.size();
    if (tp < pressureX_[0])
        FatalErrorInFunction << "Local pressure " << tp << " lower than minimum value "
        << pressureX_[0] << abort(FatalError);
    if (tp > pressureX_[length-1])
        FatalErrorInFunction << "Local pressure " << tp << " higher than maximum value "
        << pressureX_[length-1] << abort(FatalError);
    label lo = 0;
    label hi = length-1;
    while (hi - lo > 1)
    {
        label mid = lo + (hi-lo)/2;
        if (pressureX_[mid] < tp)
            lo = mid;
        else if (pressureX_[mid] > tp)
            hi = mid;
        else
            return temperatureY_[mid];
    }
    scalar lCoeff = (pressureX_[hi] - tp) / (pressureX_[hi] - pressureX_[lo]);
    scalar hCoeff = (tp - pressureX_[lo]) / (pressureX_[hi] - pressureX_[lo]);
    
    return lCoeff*temperatureY_[lo] + hCoeff*temperatureY_[hi];
}

Foam::scalar Foam::twoPhaseMixtureThermo::pv(const scalar& tp, const scalar& tT) const
{
    static const word liquidPhaseName("liquid");

    if (phase1Name() != liquidPhaseName)
        FatalErrorInFunction << "phase-1 has to be liquid" << abort(FatalError);


    static const heRhoThermopureMixtureliquidProperties& thermo =
        mesh_.lookupObject<heRhoThermopureMixtureliquidProperties>
        (
             IOobject::groupName(basicThermo::dictName, liquidPhaseName)
        );

    static const Foam::liquidProperties& liquid = thermo.mixture().properties();

    return liquid.pv(tp, tT);
}

Foam::scalar Foam::twoPhaseMixtureThermo::hl(const scalar& tp, const scalar& tT) const
{
    static const word liquidPhaseName("liquid");

    if (phase1Name() != liquidPhaseName)
        FatalErrorInFunction << "phase-1 has to be liquid" << abort(FatalError);


    static const heRhoThermopureMixtureliquidProperties& thermo =
        mesh_.lookupObject<heRhoThermopureMixtureliquidProperties>
        (
             IOobject::groupName(basicThermo::dictName, liquidPhaseName)
        );

    static const Foam::liquidProperties& liquid = thermo.mixture().properties();

    return liquid.hl(tp, tT);
}


void Foam::twoPhaseMixtureThermo::correctThermo()
{
    thermo1_->T() = T_;
    thermo1_->he() = thermo1_->he(p_, T_);
    thermo1_->correct();

    thermo2_->T() = T_;
    thermo2_->he() = thermo2_->he(p_, T_);
    thermo2_->correct();
}


void Foam::twoPhaseMixtureThermo::correct()
{
    psi_ = alpha1()*thermo1_->psi() + alpha2()*thermo2_->psi();
    mu_ = alpha1()*thermo1_->mu() + alpha2()*thermo2_->mu();
    alpha_ = alpha1()*thermo1_->alpha() + alpha2()*thermo2_->alpha();

    interfaceProperties::correct();
}


Foam::word Foam::twoPhaseMixtureThermo::thermoName() const
{
    return thermo1_->thermoName() + ',' + thermo2_->thermoName();
}


bool Foam::twoPhaseMixtureThermo::incompressible() const
{
    return thermo1_->incompressible() && thermo2_->incompressible();
}


bool Foam::twoPhaseMixtureThermo::isochoric() const
{
    return thermo1_->isochoric() && thermo2_->isochoric();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return alpha1()*thermo1_->he(p, T) + alpha2()*thermo2_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(alpha1(), cells)*thermo1_->he(T, cells)
      + scalarField(alpha2(), cells)*thermo2_->he(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::he
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->he(T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->he(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::hs() const
{
    return alpha1()*thermo1_->hs() + alpha2()*thermo2_->hs();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return alpha1()*thermo1_->hs(p, T) + alpha2()*thermo2_->hs(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(alpha1(), cells)*thermo1_->hs(T, cells)
      + scalarField(alpha2(), cells)*thermo2_->hs(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->hs(T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->hs(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::ha() const
{
    return alpha1()*thermo1_->ha() + alpha2()*thermo2_->ha();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return alpha1()*thermo1_->ha(p, T) + alpha2()*thermo2_->ha(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(alpha1(), cells)*thermo1_->ha(T, cells)
      + scalarField(alpha2(), cells)*thermo2_->ha(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::ha
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->ha(T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->ha(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::hc() const
{
    return alpha1()*thermo1_->hc() + alpha2()*thermo2_->hc();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::Cp() const
{
    return alpha1()*thermo1_->Cp() + alpha2()*thermo2_->Cp();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cp(T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cp(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::Cv() const
{
    return alpha1()*thermo1_->Cv() + alpha2()*thermo2_->Cv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cv(T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cv(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::gamma() const
{
    return alpha1()*thermo1_->gamma() + alpha2()*thermo2_->gamma();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::gamma
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->gamma(T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->gamma(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::Cpv() const
{
    return alpha1()*thermo1_->Cpv() + alpha2()*thermo2_->Cpv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cpv(T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cpv(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::CpByCpv() const
{
    return
        alpha1()*thermo1_->CpByCpv()
      + alpha2()*thermo2_->CpByCpv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::CpByCpv
(
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->CpByCpv(T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->CpByCpv(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::W() const
{
    return alpha1()*thermo1_->W() + alpha2()*thermo1_->W();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::W
(
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->W(patchi)
      + alpha2().boundaryField()[patchi]*thermo1_->W(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::nu() const
{
    return mu()/(alpha1()*thermo1_->rho() + alpha2()*thermo2_->rho());
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::nu
(
    const label patchi
) const
{
    return
        mu(patchi)
       /(
            alpha1().boundaryField()[patchi]*thermo1_->rho(patchi)
          + alpha2().boundaryField()[patchi]*thermo2_->rho(patchi)
        );
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::kappa() const
{
    return alpha1()*thermo1_->kappa() + alpha2()*thermo2_->kappa();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::kappa
(
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->kappa(patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->kappa(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::alphahe() const
{
    return
        alpha1()*thermo1_->alphahe()
      + alpha2()*thermo2_->alphahe();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::alphahe
(
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->alphahe(patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->alphahe(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    return
        alpha1()*thermo1_->kappaEff(alphat)
      + alpha2()*thermo2_->kappaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->kappaEff(alphat, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->kappaEff(alphat, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    return
        alpha1()*thermo1_->alphaEff(alphat)
      + alpha2()*thermo2_->alphaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseMixtureThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->alphaEff(alphat, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->alphaEff(alphat, patchi);
}


bool Foam::twoPhaseMixtureThermo::read()
{
    if (psiThermo::read())
    {
        return interfaceProperties::read();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
