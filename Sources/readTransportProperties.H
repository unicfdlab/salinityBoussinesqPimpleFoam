    singlePhaseTransportModel laminarTransport(U, phi);
    
    // Thermal expansion coefficient [1/K]
    dimensionedScalar beta
    (
        "beta",
        dimensionSet(0,0,0,-1,0),
        0.0,
        laminarTransport
    );
    
    // Reference temperature [K]
    dimensionedScalar SARef
    (
        "SARef",
        dimensionSet(0,0,0,1,0),
        1.0,
        laminarTransport
    );
    
    //
    dimensionedScalar DS
    (
        "DS",
        dimensionSet(0,2,-1,0,0),
        0.0,
        laminarTransport
    );
    
    //
    dimensionedScalar rho0
    (
        "rho0",
        dimensionSet(1,-3,0,0,0),
        1.0,
        laminarTransport
    );
    
    //-
    Switch salinityFVM (laminarTransport.lookup("salinityFVM"));
    
    //-
    Switch correctSalinity (laminarTransport.lookup("correctSalinity"));
    
    //-
    dimensionedScalar Sct
    (
        "Sct",
        dimless,
        1.0,
        laminarTransport
    );

    /*
    // Laminar Prandtl number
    dimensionedScalar Pr(laminarTransport.lookup("Pr"));

    // Turbulent Prandtl number
    dimensionedScalar Prt(laminarTransport.lookup("Prt"));
    */
