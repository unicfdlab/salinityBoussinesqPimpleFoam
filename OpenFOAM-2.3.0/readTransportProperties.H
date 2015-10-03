    singlePhaseTransportModel laminarTransport(U, phi);

    // Thermal expansion coefficient [1/K]
    dimensionedScalar beta(laminarTransport.lookup("beta"));

    // Reference temperature [K]
    dimensionedScalar SARef(laminarTransport.lookup("SARef"));

    //
    dimensionedScalar DS (laminarTransport.lookup("DS"));
    
    //
    dimensionedScalar rho0 (laminarTransport.lookup("rho0"));
    
    //-
    Switch salinityFVM (laminarTransport.lookup("salinityFVM"));
    
    //-
    Switch correctSalinity (laminarTransport.lookup("correctSalinity"));
    
    //-
    dimensionedScalar Sct (laminarTransport.lookup("Sct"));

    /*
    // Laminar Prandtl number
    dimensionedScalar Pr(laminarTransport.lookup("Pr"));

    // Turbulent Prandtl number
    dimensionedScalar Prt(laminarTransport.lookup("Prt"));
    */
