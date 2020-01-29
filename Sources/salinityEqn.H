{
    if (salinityFVM)
    {
        volScalarField DSEff ("DSEff", turbulence->nut()/Sct + DS);

        fvScalarMatrix salinityEqn
        (
            fvm::ddt(SA)
            + fvm::div(phi, SA)
            - fvm::laplacian(DSEff, SA)
        );

        salinityEqn.relax();
        salinityEqn.solve();
    }
    else
    {
        word salinityScheme("div(phi,SA)");

        surfaceScalarField phiSA
        (
            fvc::flux
            (
                phi,
                SA,
                salinityScheme
            )
        );

        upwind<scalar> UD (mesh, phi);
        surfaceScalarField phiSA_upwind = 
            fv::gaussConvectionScheme<scalar>(mesh,phi,UD).flux(phi,SA);
        
        MULES::explicitSolve
        (
            geometricOneField(),
            SA,
            phiSA_upwind,
            phiSA,
            zeroField(),
            zeroField(),
            oneField(),
            zeroField()
        );
    }

    dimensionedScalar domainSalinity = fvc::domainIntegrate(SA);
    Info << "Domain salinity: " << domainSalinity.value() << endl;

    rhok = 1.0 - beta*(SA - SARef);
}
