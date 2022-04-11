#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
twoPhaseMixtureThermo mixture(U, phi);
// Info<< "Mixture pv(0.0, 239.65): " << mixture.pv(0.0, 239.65) << endl;

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

Info<< "Reading thermophysical properties\n" << endl;

const volScalarField& rho1 = mixture.thermo1().rho();
const volScalarField& rho2 = mixture.thermo2().rho();

volScalarField rhoRef
(
    IOobject
    (
        "rhoRef",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);

// Phase change parameters
const scalar Cc(mixture.lookupOrDefault("Cc", 0.0));
const scalar Ce(mixture.lookupOrDefault("Ce", 0.0));
const scalar condThresh(mixture.lookupOrDefault("condThresh", 0.9));
const scalar evapThresh(mixture.lookupOrDefault("evapThresh", 0.1));
const scalar relexFac(mixture.lookupOrDefault("relaxFac", 1.0));
const scalar deltaSuperheat(mixture.lookupOrDefault("deltaSuperheat", 0.0));
Foam::Info<< "Coefficient of condensation/evaporation: " << Cc << " " << Ce << endl;

volScalarField hlField
(
    IOobject
    (
        "hlField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("hlValue", dimEnergy/dimMass, 0.0)
);

volScalarField TsatField
(
    IOobject
    (
        "TsatField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("TsatValue", dimTemperature, 0.0)
);


/*
volScalarField interfaceField
(
    IOobject
    (
        "interfaceField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("interfaceValue", dimless, 0.0)
);

volScalarField wallField
(
    IOobject
    (
        "wallField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("wallValue", dimless, 0.0)
);
*/

// mdot [kg/m3 s] mass transfer from phase2 to phase1
volScalarField mdot
(
    IOobject
    (
        "mdot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("mdot", dimMass/dimVolume/dimTime, 0.0)
);

volScalarField qvdot
(
    IOobject
    (
        "qvdot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("qv", dimEnergy/dimVolume/dimTime, 0.0)
);

volScalarField::Internal TEqnSp
(
    IOobject
    (
        "TEqnSp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimMass/dimVolume/dimTime, 0.0)
);

volScalarField::Internal TEqnSu
(
    IOobject
    (
        "TEqnSu",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimMass*dimTemperature/dimVolume/dimTime, 0.0)
);

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    mixture
);

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha1 solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField dgdt(alpha1*fvc::div(phi));

#include "createAlphaFluxes.H"

// Construct compressible turbulence model
compressibleInterPhaseTransportModel turbulence
(
    rho,
    U,
    phi,
    rhoPhi,
    alphaPhi10,
    mixture
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

#include "createMRF.H"
#include "createFvOptions.H"
