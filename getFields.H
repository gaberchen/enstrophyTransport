    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject UMeanheader
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );


    if (Uheader.headerOk())
    {
        Info<< "    Reading U" << endl;
	volVectorField U(Uheader, mesh);

	Info<< "    Reading UMean" << endl;
        volVectorField UMean(UMeanheader, mesh);

	volVectorField UPrime(U-UMean);
	volVectorField vorticityPrime(fvc::curl(UPrime));
        volVectorField vorticityMean(fvc::curl(UMean));
	volTensorField UMeanGrad(fvc::grad(UMean));
	volTensorField UPrimeGrad(fvc::grad(UPrime));
	volTensorField vorticityMeanGrad(fvc::grad(vorticityMean));
	volTensorField vorticityPrimeGrad(fvc::grad(vorticityPrime));
	
	Info<< "    Calculating P1 term" << endl;
        volScalarField enstrophyProd1
	(
	    IOobject
	    (
		"enstrophyProd1",
	    	runTime.timeName(),
	    	mesh,
	    	IOobject::NO_READ
	    ),
	    (vorticityPrime * vorticityPrime) && UMeanGrad
	);
	
	Info<< "    Calculating P2 term" << endl;
        volScalarField enstrophyProd2
	(
	    IOobject
	    (
		"enstrophyProd2",
	    	runTime.timeName(),
	    	mesh,
	    	IOobject::NO_READ
	    ),
	    (vorticityPrime & UPrimeGrad) & vorticityMean
	);

	Info<< "    Calculating P3 term" << endl;
	volScalarField enstrophyProd3
	(
	    IOobject
	    (
		"enstrophyProd3",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ
	    ),
	    (UPrime * vorticityPrime) && vorticityMeanGrad
	);

	Info<< "    Calculating P4 term" << endl;
	volScalarField enstrophyProd4
	(
	    IOobject
	    (
		"enstrophyProd4",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ
	    ),
	    (vorticityPrime * vorticityPrime) && UPrimeGrad
	); 

	Info<< "    Calculating transport term" << endl;
	volScalarField enstrophyTransport
	(
	    IOobject
	    (
		"enstrophyTransort",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ
	    ),
	    -0.5 * fvc::div( UPrime & (vorticityPrime * vorticityPrime))
	);

        IOdictionary transportProperties
        (
             IOobject
             (
                 "transportProperties",
                 runTime.constant(),
                 mesh,
                 IOobject::MUST_READ_IF_MODIFIED,
                 IOobject::NO_WRITE
             )
        );

        dimensionedScalar nu(transportProperties.lookup("nu"));

	Info<< "    Calculating viscous diffusion term" << endl;
	volScalarField enstrophyDiffusion
	(
	    IOobject
	    (
		"enstrophyDiffusion",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ
	    ),
	    nu * fvc::laplacian( 0.5*magSqr(vorticityPrime) )
	);

	Info<< "    Calculating viscous dissipation term" << endl;
	volScalarField enstrophyDissipation
	(
	    IOobject
	    (
		"enstrophyDissipation",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ
	    ),
	    -1.0 * nu * (vorticityPrimeGrad && vorticityPrimeGrad)
	);

    if (writeResults)
    {
	Info<< "    Writing terms to files" << endl;        
	enstrophyProd1.write();
        enstrophyProd2.write();
	enstrophyProd3.write();
	enstrophyProd4.write();
	enstrophyTransport.write();
	enstrophyDiffusion.write();
	enstrophyDissipation.write();
    }

    }
    else
    {
        Info<< "    No UMean" << endl;
    }

