{
    if (salinityFVM)
    {
	/*
	kappat = turbulence->nut()/Prt;
	kappat.correctBoundaryConditions();

	volScalarField kappaEff("kappaEff", turbulence->nu()/Pr + kappat);
	*/
	
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
	word alpharScheme("div(phirb,alpha)");

        surfaceScalarField phiSA
        (
            fvc::flux
            (
                phi,
                SA,
                salinityScheme
            )
//don't use compressive flux
//            + fvc::flux
//            (
//                -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
//                alpha1,
//                alpharScheme
//            )
        );
        
        volScalarField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            fvc::div(phi) * 0.0
        );

        volScalarField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
	    -Sp*SA
        );
        
        MULES::implicitSolve
        (
            geometricOneField(),
            SA,
            phi,
            phiSA,
            Sp,
            Su,
            1,
            0
        );
    }

    dimensionedScalar domainSalinity = fvc::domainIntegrate(SA);
    Info << "Domain salinity: " << domainSalinity.value() << endl;

    rhok = 1.0 - beta*(SA - SARef);
}
