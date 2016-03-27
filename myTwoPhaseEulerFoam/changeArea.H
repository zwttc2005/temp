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
    	if (celli >= (celln-50) && celli < (celln))
        {
		        celld = (double) (celli-(celln-100))/30;  
		        scaleFactor = Foam::tanh(celld);	
		        flowArea[celli] = flowArea[(celln-100)]-flowArea[(celln-100)]*scalar(scaleFactor);	
        }  
        /*else if (celli >= (celln-4) && celli < (celln))
	{
		        celld = (double) (celli-(celln-52))/200;  
		        scaleFactor = Foam::tanh(celld);	
		        flowArea[celli] = flowArea[(celln-52)]+flowArea[(celln-52)]*scalar(scaleFactor);	
	}*/
	else if (celli >= celln)
	{
		        flowArea[celli] = flowArea[(celln-1)];
	}
    	
    }

    Info << "flowArea[end] = " << flowArea[celln] << nl << endl;
    Info << "flowArea[nozzle] = " << flowArea[(celln-52)] << nl << endl;



    // assign flowAreaGrad to each cell centre
    forAll(flowAreaGrad, celli)
    {
	   if  (celli < (celln-100))    
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


