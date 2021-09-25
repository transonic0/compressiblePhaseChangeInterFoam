volScalarField::Internal Sp
(
    IOobject
    (
        "Sp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dgdt.dimensions(), 0)
);

volScalarField::Internal Su
(
    IOobject
    (
        "Su",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dgdt.dimensions(), 0)
);

// Phase change mass transfer
const volScalarField::Internal& tmpPCells = p.internalField();
const volScalarField::Internal& tmpTCells = T.internalField();
volScalarField::Internal& mdotCells = mdot.ref();
volScalarField::Internal& qvdotCells = qvdot.ref();
volScalarField interfacialArea = mag( fvc::grad(alpha1) );
// const volScalarField::Internal& tmpArea = interfacialArea.internalField();
dimensionedScalar totalArea = fvc::domainIntegrate(interfacialArea);
Foam::Info<< "System total interfacial area (m2):  " << totalArea.value() << endl;
// volScalarField mixtureCp = mixture.Cp();
// const volScalarField::Internal& mixtureCpCells = mixture.Cp()().internalField();
const volScalarField mixtureCp =
    alpha1 * mixture.thermo1().Cp()() + alpha2 * mixture.thermo2().Cp()();
const volScalarField::Internal& mixtureCpCells = mixtureCp.internalField();

// for (label i = 0; i < 1001; ++i)
// {
//     const scalar p0 = 1.0 + 2e07/1000*i;
//     Foam::Info<< p0 << "," << mixture.Tsat(p0) << endl;
// }

forAll(mdotCells, celli)
{
    const scalar pc = tmpPCells[celli];
    const scalar Tc = tmpTCells[celli];
    const scalar Tsat = mixture.Tsat(pc);
    const scalar hl = mixture.hl(pc, Tc);
    const scalar cpc = mixtureCpCells[celli];
    // Foam::Info << "Tc Tsat cp hl: " << Tc << " " << Tsat << " " << cpc << " " << hl << endl;
    // const scalar psat = mixture.pv(pc, Tc);
    // const scalar tmpCoeff = Foam::sqrt( MWT / (2.0*3.14159265359*8.314462618153*Tc));
    // // Foam::Info<< "pc Tc psat tmpCoeff A: " << pc << " " << Tc << " " << psat << " "
    // //     << tmpCoeff << " " << interfacialArea.internalField()[celli] << endl;
    // mdotCells[celli] = tmpCoeff * (Cc*pc - Ce*psat) * tmpArea[celli];
    // Su[celli] += mdotCells[celli] * (alpha2[celli]/rho1[celli] + alpha1[celli]/rho2[celli]);

    scalar Tr = 0.0;
    if (Tc > Tsat)  // Evaporation
    {
        Tr = (Tc - Tsat) / Tsat;
        mdotCells[celli] = -Ce*alpha1[celli]*rho1[celli]*Tr;
        Sp[celli] -= Ce * Tr * (alpha2[celli] + rho1[celli]/rho2[celli] * alpha1[celli]);

        qvdotCells[celli] = -Ce*alpha1[celli]*rho1[celli]*hl*Tr;
        TEqnSp[celli] = -Ce*alpha1[celli]*rho1[celli]*hl/(Tsat*cpc);
        TEqnSu[celli] = Ce*alpha1[celli]*rho1[celli]*hl/(cpc);
    }
    else if (Tc < Tsat)  // Condensation
    {
        Tr = (Tsat - Tc) / Tsat;
        mdotCells[celli] = Cc*alpha2[celli]*rho2[celli]*Tr;
        Sp[celli] -= Cc * Tr * (alpha1[celli] + rho2[celli]/rho1[celli] * alpha2[celli]);
        Su[celli] += Cc * Tr * (alpha1[celli] + rho2[celli]/rho1[celli] * alpha2[celli]);

        qvdotCells[celli] = Cc*alpha2[celli]*rho2[celli]*hl*Tr;
        TEqnSp[celli] = -Cc*alpha2[celli]*rho2[celli]*hl/(Tsat*cpc);
        TEqnSu[celli] = Cc*alpha2[celli]*rho2[celli]*hl/(cpc);
    }
    else  // No phase change
    {
        mdotCells[celli] = 0.0;
        qvdotCells[celli] = 0.0;
        TEqnSp[celli] = 0.0;
        TEqnSu[celli] = 0.0;
    }
}

Foam::Info<< "mdot min/max: " << min(mdot).value() << " " << max(mdot).value() << endl;
Foam::Info<< "qvdot min/max: " << min(qvdot).value() << " " << max(qvdot).value() << endl;

forAll(dgdt, celli)
{
    if (dgdt[celli] > 0.0)
    {
        Sp[celli] -= dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
        Su[celli] += dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
    }
    else if (dgdt[celli] < 0.0)
    {
        Sp[celli] += dgdt[celli]/max(alpha1[celli], 1e-4);
    }
}

volScalarField::Internal divU
(
    mesh.moving()
  ? fvc::div(phiCN() + mesh.phi())
  : fvc::div(phiCN())
);
