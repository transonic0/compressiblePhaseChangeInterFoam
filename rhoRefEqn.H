{
    fvScalarMatrix rhoRefEqn
    (
        fvm::ddt(rhoRef)
      + fvc::div(rhoPhi)
    );

    rhoRefEqn.solve();

    Info<< "rhoRef min/max: " << min(rhoRef).value()
        << " " << max(rhoRef).value() << endl;
    Info<< "p min/max: " << min(p).value()
        << " " << max(p).value() << endl;
}
