/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Global
    continuityErrs

Description
    Calculates and prints the continuity errors.

\*---------------------------------------------------------------------------*/

{
    #include "rhoRefEqn.H"

    dimensionedScalar totalMass = fvc::domainIntegrate(rhoRef);
    Foam::Info<< "System total mass (kg):  " << totalMass.value() << endl;
    static scalar mass0 = totalMass.value();
    Foam::Info<< "System total mass deviation (%):  " << (totalMass.value() - mass0)/mass0 << endl;

    dimensionedScalar totalPhaseChange = fvc::domainIntegrate(mdot);
    Foam::Info<< "System total phase change rate (kg/s):  " << totalPhaseChange.value() << endl;

    // rho = alpha1*rho1 + alpha2*rho2;

    scalar sumLocalContErr =
        (fvc::domainIntegrate(mag(rhoRef - rho))/totalMass).value();

    scalar globalContErr =
        (fvc::domainIntegrate(rhoRef - rho)/totalMass).value();

    cumulativeContErr += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr
        << endl;
}

// ************************************************************************* //
