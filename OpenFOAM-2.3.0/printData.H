{
    volScalarField rho = rhok*rho0;
    gradRhoTilde = fvc::grad(rho - initialRho);
    rho.rename("rho");
    if (runTime.outputTime())
    {
	rho.write();
	gradRhoTilde.write();
    }
}

