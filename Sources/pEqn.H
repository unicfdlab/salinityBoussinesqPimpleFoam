{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rAU));
    
    U = HbyA;
    
    bool ddtPhiCorr = false;
    if (pimple.dict().found("ddtPhiCorr"))
    {
        pimple.dict().lookup("ddtPhiCorr") >> ddtPhiCorr;
    }
    
    phi = (fvc::interpolate(U) & mesh.Sf());
    
    if (ddtPhiCorr)
	phi += rAUf * fvc::ddtCorr(U, phi);

    //surfaceScalarField buoyancyPhi(rAUf*ghf*fvc::snGrad(rhokstar)*mesh.magSf());
    surfaceScalarField buoyancyPhi(-rAUf*(linearInterpolate(g*rhokstar)&mesh.Sf()));
    phi -= buoyancyPhi;
    
    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phi, rAUf);
    
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phi)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi -= p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            //U -= rAU*fvc::reconstruct((buoyancyPhi + p_rghEqn.flux())/rAUf);
            U -= (fvc::grad(p_rgh) - rhokstar*g)*rAU;
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"
    
    if (correctSalinity)
    {
	#include "salinityEqn.H"
    }

    //p = p_rgh + rhok*gh;
    //if (p_rgh.needReference())
    //{
    //    p += dimensionedScalar
    //    (
    //        "p",
    //        p.dimensions(),
    //        pRefValue - getRefCellValue(p, pRefCell)
    //    );
    //    p_rgh = p - rhok*gh;
    //}
}
