    //constant and constant vectors
		    
    double celld; 
    
    double scaleFactor;

    double xDim;	

    int celln;

    // declare cell-to-cell length
    volScalarField xDimDim // cell-to-cell length	
    (
     IOobject
     (

      "xDimDim",
      runTime.timeName(),
      mesh

     ),

     mesh,
     dimensionedScalar("0",dimensionSet (0,1,0,0,0),0)
     
    );	
    // declare the cell centre variable 
    volScalarField flowArea // flow area of that specific cell 	
    (
     IOobject
     (

      "flowArea",
      runTime.timeName(),
      mesh

     ),

     mesh,
     dimensionedScalar("0",dimensionSet (0,2,0,0,0),0)
     
    );	
   

    // declare the gradient of area along flow direction
    volScalarField flowAreaGrad 	
    (
     IOobject
     (

      "flowAreaGrad",
      runTime.timeName(),
      mesh

     ),

     mesh,
     dimensionedScalar("0",dimensionSet (0,1,0,0,0),0)
     
    );	
    
    // declare the area part of the source term
    volScalarField areaSource 	
    (
     IOobject
     (

      "areaSource",
      runTime.timeName(),
      mesh

     ),

     mesh,
     dimensionedScalar("0",dimensionSet (0,-1,0,0,0),0)
     
    );	
    
    // declare bulk density 
    volScalarField rho_bulk 	
    (
     IOobject
     (

      "rho_bulk",
      runTime.timeName(),
      mesh

     ),

     mesh,
     dimensionedScalar("0",dimensionSet (1,-3,0,0,0),0)
     
    );	
    
    // declare bulk velocity 
    volScalarField U_bulk 	
    (
     IOobject
     (

      "U_bulk",
      runTime.timeName(),
      mesh

     ),

     mesh,
     dimensionedScalar("0",dimensionSet (0,1,-1,0,0),0)
     
    );	
   
    // declare bulk enthalpy 
    volScalarField he_bulk 	
    (
     IOobject
     (

      "he_bulk",
      runTime.timeName(),
      mesh

     ),

     mesh,
     dimensionedScalar("0",dimensionSet (0,2,-2,0,0),0)
     
    );	

    // declare volume weighted compressibility 
    volScalarField psi_bulk 	
    (
     IOobject
     (

      "psi_bulk",
      runTime.timeName(),
      mesh

     ),

     mesh,
     dimensionedScalar("0",dimensionSet (0,-2,2,0,0),0)
     
    );	

     // obtain cell length 
    const faceList & ff = mesh.faces();
    const pointField & pp = mesh.points();

    forAll(mesh.C(), celli)
    {
    	const cell & cc = mesh.cells()[celli];
	labelList pLabels(cc.labels(ff));
	pointField pLocal(pLabels.size(), vector::zero);

	forAll (pLabels, pointi)
	{
		pLocal[pointi] = pp[pLabels[pointi]];
	}

	xDim = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
	xDimDim[celli] = scalar(xDim);
    }


    // assign flowarea value to the cell centre 
    // obtain number of cells as well
    celln = -1;
    forAll(flowArea,celli)
    {
    	flowArea[celli] = 1; // initialisation
    	celln += 1;
    }

    forAll(flowArea,celli)
    {
    	if (celli >= (celln-100) && celli < (celln))
	{
		celld = (double) (celli-(celln-100))/50;  
		scaleFactor = Foam::tanh(celld);	
		flowArea[celli] = flowArea[(celln-100)]-flowArea[(celln-100)]*scalar(scaleFactor);	
	}
	//else if (celli > (celln-50) && celli < celln)
	//{
		//celld = (double) (celli-(celln-50))/100;
		//scaleFactor = Foam::tanh(celld);
		//flowArea[celli] =flowArea[(celln-50)]+5.0*scalar(scaleFactor);//flowArea[(celln-50)]*scalar(scaleFactor);
	//}
	else if (celli >= celln)
	{
		flowArea[celli] = flowArea[(celln-1)];
	}
    	
    }

    // assign flowAreaGrad to each cell centre
    forAll(flowAreaGrad, celli)
    {
	if (celli < (celln-100))    
	{
		flowAreaGrad[celli] = 0; 
    	}
    	else if (celli >= (celln-100) && celli < celln) 
	{
		flowAreaGrad[celli] = (flowArea[celli+1] - flowArea[celli])/xDimDim[celli];
	}
	else if (celli >= celln)
	{	
		flowAreaGrad[celli] = 0;
	}	
    }
	
    // assign areaSource to each cell centre
    forAll(areaSource, celli)
    {
	if (celli < (celln-100))    
	{
		areaSource[celli] = 0; 
    	}
    	else if (celli >= (celln-100) && celli < celln) 
	{
		areaSource[celli] = flowAreaGrad[celli]/flowArea[celli];
	}
	else if (celli >= celln)
	{	
		areaSource[celli] = 0;
	}	
    }


