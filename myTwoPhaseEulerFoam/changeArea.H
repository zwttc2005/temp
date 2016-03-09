    //constant and constant vectors
		    
    double celld; 
    
    double scaleFactor;

    double xDim;

    double xDimSum;

    int celln;
    
    // obtain number of cells
    // initialise flowArea 
    celln = -1;
    forAll(flowArea,celli)
    {
    	flowArea[celli] = 1; // initialisation
    	celln += 1;
    }
	
    // obtain cell length 
    const faceList & ff = mesh.faces();
    const pointField & pp = mesh.points();
    static int flagLength = 0;	

    xDimSum = 0.0;	

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
    	
	if (celli <= celln) xDimSum = xDimSum+xDim;
	if (xDimSum >= 256.0 && flagLength == 0)
	{
		Info << "Pipe end cell index = " << celli << nl << endl;
		flagLength++;
	}
    }

    //Info << "Dimension = " << xDimSum << nl << endl; 	


    // assign flowarea value to the cell centre 


    forAll(flowArea,celli)
    {
    	if (celli >= (celln-51) && celli < (celln))
	{
		celld = (double) (celli-(celln-51))/200;  
		scaleFactor = Foam::tanh(celld);	
		flowArea[celli] = flowArea[(celln-51)]+flowArea[(celln-51)]*scalar(scaleFactor);	
	}
	//else if (celli > (celln-50) && celli < celln)
	//{
		//celld = (double) (celli-(celln-50))/100;
		//scaleFactor = Foam::tanh(celld);
		//flowArea[celli] =flowArea[(celln-50)]+5.0*scalar(scaleFactor);//flowArea[(celln-50)]*scalar(scaleFactor);
	//}
	else if (celli >= celln-51)
	{
		flowArea[celli] = flowArea[(celln-51)];
	}
    	
    }

    // assign flowAreaGrad to each cell centre
    forAll(flowAreaGrad, celli)
    {
	if (celli < (celln-51))    
	{
		flowAreaGrad[celli] = 0; 
    	}
    	else if (celli >= (celln-51) && celli < celln) 
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
	if (celli < (celln-51))    
	{
		areaSource[celli] = 0; 
    	}
    	else if (celli >= (celln-51) && celli < celln) 
	{
		areaSource[celli] = flowAreaGrad[celli]/flowArea[celli];
	}
	else if (celli >= celln-51)
	{	
		areaSource[celli] = 0;
	}	
    }


