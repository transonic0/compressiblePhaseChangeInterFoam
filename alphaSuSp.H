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

// calculate Tsat and hl
scalarField& TsatFieldCells = TsatField.primitiveFieldRef();
scalarField& hlFieldCells = hlField.primitiveFieldRef();
forAll(tmpPCells, celli)
{
    TsatFieldCells[celli] = mixture.Tsat( tmpPCells[celli] );
    hlFieldCells[celli] = mixture.hl( tmpPCells[celli], tmpTCells[celli] );
}
forAll(p.boundaryField(), patchI)
{
    const fvPatchScalarField& pp = p.boundaryField()[patchI];
    const fvPatchScalarField& pT = T.boundaryField()[patchI];
    fvPatchScalarField& pTsat = TsatField.boundaryFieldRef()[patchI];
    fvPatchScalarField& phl = hlField.boundaryFieldRef()[patchI];
    forAll(pp, facei)
    {
        pTsat[facei] = mixture.Tsat( pp[facei] );
        phl[facei] = mixture.hl( pp[facei], pT[facei] );
    }
}

const volScalarField heField = mixture.he(p, T);
// const volScalarField heSatField = mixture.he(p, TsatField);

// for (label i = 0; i < 1001; ++i)
// {
//     const scalar p0 = 1.0 + 2e07/1000*i;
//     Foam::Info<< p0 << "," << mixture.Tsat(p0) << endl;
// }


/**
 * 20220323 new interface equilibrium phase change model
 *
 */

/*
std::vector<MeshGraph::CellFacePair> condIntCellFacePairs;
std::vector<MeshGraph::CellFacePair> evapIntCellFacePairs;

interfaceMeshGraph.Reset();
interfaceMeshGraph.GetInterfaceCellFacePairs(condIntCellFacePairs, condThresh);
interfaceMeshGraph.Reset();
interfaceMeshGraph.GetInterfaceCellFacePairs(evapIntCellFacePairs, evapThresh);

surfaceScalarField Tface = fvc::interpolate(T);
surfaceScalarField pface = fvc::interpolate(p);
surfaceScalarField TsatFace = fvc::interpolate(TsatField);
interfaceField = 0.0;

for
(
    std::vector<MeshGraph::CellFacePair>::iterator it = condIntCellFacePairs.begin();
    it != condIntCellFacePairs.end();
    it++
)
{
    if (Tface[(*it).f] <= TsatFace[(*it).f])
    {
        interfaceField[ (*it).c1 ] = 1;
        interfaceField[ (*it).c2 ] = 1;
    }
}

for
(
    std::vector<MeshGraph::CellFacePair>::iterator it = evapIntCellFacePairs.begin();
    it != evapIntCellFacePairs.end();
    it++
)
{
    if (Tface[(*it).f] >= TsatFace[(*it).f])
    {
        interfaceField[ (*it).c1 ] = 1;
        interfaceField[ (*it).c2 ] = 1;
    }
}

labelList wallCells;
forAll (mesh.boundary(), pI)
{
    if ( isA<wallFvPatch>(mesh.boundary()[pI]) )
        wallCells.append( mesh.boundary()[pI].faceCells() );
}
wallField = 0.0;
forAll ( wallCells, cI )
{
    wallField[ wallCells[cI] ] = 1;
    interfaceField[ wallCells[cI] ] = 1;
}

mdot = dimensionedScalar( "dummy", dimMass/dimVolume/dimTime, 0.0 );
qvdot = dimensionedScalar( "dummy", dimEnergy/dimVolume/dimTime, 0.0 );
TEqnSp = dimensionedScalar("dummy", dimMass/dimVolume/dimTime, 0.0 );
TEqnSu = dimensionedScalar( "dummy", dimMass*dimTemperature/dimVolume/dimTime, 0.0);

const dimensionedScalar dT = runTime.deltaT();

qvdot = interfaceField * rho * mixtureCp * (T - TsatField) / dT;
Foam::Info<< "* Unlimited * qvdot min/max: " << min(qvdot).value() << " " << max(qvdot).value() << endl;
qvdot = interfaceField * rho * (heField - heSatField) / dT;
Foam::Info<< "**Unlimited** qvdot min/max: " << min(qvdot).value() << " " << max(qvdot).value() << endl;
volScalarField limEvap = (1.0 - wallField) * alpha1 * rho1 * hlField / dT;
volScalarField limCond = (1.0 - alpha1) * rho2 * hlField / dT;
volScalarField qvdot_fluid = neg(qvdot) * max(qvdot, -limCond) + pos(qvdot) * min(qvdot, limEvap);
volScalarField PCV_fac = dT * (qvdot / hlField) * (scalar(1.0)/rho2 - scalar(1.0)/rho1);
volScalarField qvdot_vol = qvdot * mag( min( max( 1.0/(PCV_fac+SMALL), -1.0 ), (1.0-wallField) ) );
qvdot = pos(qvdot) * min( min(qvdot, qvdot_fluid), qvdot_vol ) + neg(qvdot) * max( max( qvdot, qvdot_fluid ), qvdot_vol );
qvdot = -relexFac * qvdot;
mdot = qvdot / hlField;

forAll(mdotCells, celli)
{
    Su[celli] = (alpha1[celli]/rho2[celli] + alpha2[celli]/rho1[celli]) * mdot[celli];

    TEqnSu[celli] = qvdot[celli] / mixtureCpCells[celli];
}
Foam::Info<< "qvdot min/max: " << min(qvdot).value() << " " << max(qvdot).value() << endl;
Foam::Info<< "mdot min/max: " << min(mdot).value() << " " << max(mdot).value() << endl;
*/


/**
 * Original Lee phase change model
 *
 */

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
        /**
         * 20220411
         * Superheat Lee model
         */
        if (Tc >= Tsat + deltaSuperheat)
        {
            Tr = (Tc - Tsat) / Tsat;
            mdotCells[celli] = -Ce*alpha1[celli] * rho1[celli]*Tr;
            Sp[celli] -= Ce * Tr * (alpha2[celli] + rho1[celli]/rho2[celli] * alpha1[celli]);

            qvdotCells[celli] = -Ce*alpha1[celli]*rho1[celli]*hl*Tr;
            TEqnSp[celli] = -Ce*alpha1[celli]*rho1[celli]*hl/(Tsat*cpc);
            TEqnSu[celli] = Ce*alpha1[celli]*rho1[celli]*hl/(cpc);
        }
        else if (alpha1 < 1)
        {
            Tr = (Tc - Tsat) / Tsat;
            mdotCells[celli] = -Ce*alpha1[celli] * rho1[celli]*Tr;
            Sp[celli] -= Ce * Tr * (alpha2[celli] + rho1[celli]/rho2[celli] * alpha1[celli]);

            qvdotCells[celli] = -Ce*alpha1[celli]*rho1[celli]*hl*Tr;
            TEqnSp[celli] = -Ce*alpha1[celli]*rho1[celli]*hl/(Tsat*cpc);
            TEqnSu[celli] = Ce*alpha1[celli]*rho1[celli]*hl/(cpc);
        }
        else
        {
            mdotCells[celli] = 0.0;
            qvdotCells[celli] = 0.0;
            TEqnSp[celli] = 0.0;
            TEqnSu[celli] = 0.0;
        }

        /**
         * Original
         */
        /*
        Tr = (Tc - Tsat) / Tsat;
        // // 20220307
        // mdotCells[celli] = -Ce*alpha1[celli] * alpha2[celli] *rho1[celli]*Tr;
        // Sp[celli] -= Ce * Tr * (alpha2[celli] + rho1[celli]/rho2[celli] * alpha1[celli]) * alpha2[celli];

        // qvdotCells[celli] = -Ce*alpha1[celli]* alpha2[celli] *rho1[celli]*hl*Tr;
        // TEqnSp[celli] = -Ce*alpha1[celli]* alpha2[celli] *rho1[celli]*hl/(Tsat*cpc);
        // TEqnSu[celli] = Ce*alpha1[celli]* alpha2[celli] *rho1[celli]*hl/(cpc);


        //original
        mdotCells[celli] = -Ce*alpha1[celli] * rho1[celli]*Tr;
        Sp[celli] -= Ce * Tr * (alpha2[celli] + rho1[celli]/rho2[celli] * alpha1[celli]);

        qvdotCells[celli] = -Ce*alpha1[celli]*rho1[celli]*hl*Tr;
        TEqnSp[celli] = -Ce*alpha1[celli]*rho1[celli]*hl/(Tsat*cpc);
        TEqnSu[celli] = Ce*alpha1[celli]*rho1[celli]*hl/(cpc);
        */
    }
    else if (Tc < Tsat)  // Condensation
    {
        Tr = (Tsat - Tc) / Tsat;

        // //20220307
        // mdotCells[celli] = Cc*alpha1[celli]*alpha2[celli]*rho2[celli]*Tr;
        // Sp[celli] -= Cc * Tr * alpha1[celli]* (alpha1[celli] + rho2[celli]/rho1[celli] * alpha2[celli]);
        // Su[celli] += Cc * Tr * alpha1[celli]* (alpha1[celli] + rho2[celli]/rho1[celli] * alpha2[celli]);

        // qvdotCells[celli] = Cc*alpha1[celli]*alpha2[celli]*rho2[celli]*hl*Tr;
        // TEqnSp[celli] = -Cc*alpha1[celli]*alpha2[celli]*rho2[celli]*hl/(Tsat*cpc);
        // TEqnSu[celli] = Cc*alpha1[celli]*alpha2[celli]*rho2[celli]*hl/(cpc);


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

/**
 * End original Lee model
 */

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
