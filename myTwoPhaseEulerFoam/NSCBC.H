	//-----------------------------------------------------------
	// sonic/supersonic outflow boundary condition implementation
	// 
	// In 0/p file, sideRight boundary condition must be of 
	// fixed value type 
	//-----------------------------------------------------------		
	label patchID = mesh.boundaryMesh().findPatchID("sideRight");
	
	//const scalarField& p_test = p.boundaryField()[patchID];
	    
	// make a reference to cell adjacent to boundary patch 
  	scalar p_interior = p[(celln)];
	scalar U_interior = mag(U1[(celln)])*alpha1[celln]+mag(U2[celln])*alpha2[celln];
 	scalar rho_interior = alpha1[(celln)]*rho1[(celln)]+alpha2[(celln)]*rho2[(celln)];

	Info << "pressure1 = " << p_interior << nl << endl;	
	Info << "velocity1 = " << U_interior << nl << endl;
	Info << "density1 = " << rho_interior << nl << endl;	
 	Info << "temperatue1 = " << thermo1.T()[celln] << nl << endl;
	//Info << "enthalpy1 = " << thermo1.he()[celln] << nl << endl; 

  	// declare the properties which is calculated based on characteristic
  	// analysis
  	scalar p_ghost_update, U_ghost_update, rho_ghost_update, T_ghost_update, rho1_ghost_update, rho2_ghost_update, alpha1_ghost_update;
  	scalar U_ghost, p_ghost, rho_ghost, rho1_ghost, rho2_ghost, T1_ghost, T2_ghost, alpha1_ghost, alpha2_ghost;
	scalar sos, sos1, sos2;

  	// obtaining the ghost cell information

	if (runTime.time().value() == 0.0)
	{	
	    U_ghost = scalar(0);
	    p_ghost = scalar(599000);
	    rho_ghost = scalar(33.295);
      	    rho1_ghost = scalar(33.295);
      	    rho2_ghost = scalar(1088.9);
      	    T1_ghost = scalar(220.00);
      	    T2_ghost = scalar(220.00);
            alpha1_ghost = scalar(1.0)*(rho_ghost/rho1_ghost);
            alpha2_ghost = scalar(1.0) - alpha1_ghost;
	}
  	else
  	{
  	        U_ghost = U_ghost_update;       //mag(U.boundaryField()[patchID][0]);
      	    p_ghost = p_ghost_update;       //p.boundaryField()[patchID][0];
	    rho_ghost = alpha1.boundaryField()[patchID][0]*rho1.boundaryField()[patchID][0]+alpha2.boundaryField()[patchID][0]*rho2.boundaryField()[patchID][0];
	    rho1_ghost = rho1.boundaryField()[patchID][0];
      	    rho2_ghost = rho2.boundaryField()[patchID][0];
      	    T1_ghost = thermo1.T().boundaryField()[patchID][0]; 
      	    T2_ghost = thermo2.T().boundaryField()[patchID][0];
      	    alpha1_ghost = alpha1.boundaryField()[patchID][0]; 
      	    alpha2_ghost = alpha2.boundaryField()[patchID][0];
  	}
	
    	  Info << "pressure2 = " << p.boundaryField()[patchID][0] << nl << endl;
	  Info << "velocity2 = " << U1.boundaryField()[patchID][0] << nl << endl;
	  Info << "velocity2 = " << U2.boundaryField()[patchID][0] << nl << endl;
	  Info << "densityVap2 = " << rho1.boundaryField()[patchID][0] << nl << endl;
	  Info << "densityLiq2 = " << rho2.boundaryField()[patchID][0] << nl << endl;	
          Info << "temperature2 = " << thermo2.T().boundaryField()[patchID][0] << nl << endl;			


    //
    // obtaining the speed of sound      
    //
    // Thermodynamic calculations
	 
    long ierr = 0; 
    char herr[255];
    double x_cal[] = {1.0}, xliq_cal[] = {0.0}, xvap_cal[] = {0.0}, q_cal = NULL,T_cal = NULL, rho_cal = NULL, rholiq_cal=NULL, rhovap_cal=NULL, p_cal=NULL, h_cal = NULL, s_cal = NULL, e_cal = NULL, cp_cal = NULL, cv_cal = NULL, w_cal = NULL;
    double hjt = NULL;	
    double MW = 44.01; // kg/kmol

    T_cal = T1_ghost, rho_cal = rho1_ghost/MW;
    THERMdll(T_cal,rho_cal,x_cal,p_cal,e_cal,h_cal,s_cal,cv_cal,cp_cal,w_cal,hjt);
	   	
    if (w_cal = w_cal) sos1 = w_cal;
    else if (w_cal != w_cal) sos1 = sos1;	

    T_cal = T2_ghost, rho_cal = rho2_ghost/MW;
    THERMdll(T_cal,rho_cal,x_cal,p_cal,e_cal,h_cal,s_cal,cv_cal,cp_cal,w_cal,hjt);

    if (w_cal = w_cal) sos2 = w_cal;
    else if (w_cal != w_cal) sos2 = sos2;	

    // bulk speed of sound (weighted against volume fractions)
    if (alpha1_ghost > 0.0 && alpha1_ghost < 1.0)
    {
    	sos = Foam::sqrt(1.0/(rho_ghost*(alpha1_ghost/(rho1_ghost*sos1*sos1)+alpha2_ghost/(rho2_ghost*sos2*sos2)))); 
    }
    else if (alpha1_ghost == 1.0) sos = sos1;
    else if (alpha1_ghost == 0.0) sos = sos2;
  	
    Info << "sound speed =" << sos1 << nl << endl;	
    Info << "sound speed =" << sos2 << nl << endl;	
    Info << "sound speed =" << sos << nl << endl;  
    //----------------------------------------------------------			
	
  
    // obtaining current timestep deltaT
    scalar dt = runTime.deltaT().value();
		
    // wave amplitude calculation
    if (U_ghost > sos)
    {
	  	scalar wave1 = ((U_ghost-sos)/xDimDim[celln])*(p_ghost-p_interior-rho_ghost*sos*(U_ghost-U_interior));
	  	scalar wave2 = (U_ghost/xDimDim[celln])*(p_ghost - p_interior-sos*sos*(rho_ghost-rho_interior));
	  	scalar wave3 = ((U_ghost+sos)/xDimDim[celln])*(p_ghost-p_interior+rho_ghost*sos*(U_ghost-U_interior));

		Info << "wave1 = " << wave1 << nl << endl;
		Info << "wave2 = " << wave2 << nl << endl;	  
		Info << "wave3 = " << wave3 << nl << endl;

      		scalar DpDt = -0.5*(wave3+wave1);
	  	scalar DuDt = -0.5/(rho_ghost*sos)*(wave3-wave1);
	  	scalar DrhoDt = 1.0/(sos*sos)*(DpDt+wave2);

	  	p_ghost_update = p_ghost + DpDt*dt;
	  	U_ghost_update = U_ghost + DuDt*dt;
	  	rho_ghost_update = rho_ghost + DrhoDt*dt;
    }
    else if (U_ghost <= sos)
    {
	  	scalar wave2 = 0.0; //(U_ghost/xDimDim[celln])*(p_ghost - p_interior-sos*sos*(rho_ghost-rho_interior));
	  	scalar wave3 = ((U_ghost+sos)/xDimDim[celln])*(p_ghost-p_interior+rho_ghost*sos*(U_ghost-U_interior));
    		scalar wave1 = -wave3;
  
		Info << "wave1 = " << wave1 << nl << endl;
		Info << "wave2 = " << wave2 << nl << endl;	  
		Info << "wave3 = " << wave3 << nl << endl;  
	  	
		scalar DpDt = -0.5*(wave3+wave1);
	  	scalar DuDt = -0.5/(rho_ghost*sos)*(wave3-wave1);
	  	scalar DrhoDt = 1.0/(sos*sos)*(DpDt+wave2);

  		p_ghost_update = scalar(599000); //p_ghost + DpDt*dt;
	  	U_ghost_update = U_ghost + DuDt*dt;
	  	rho_ghost_update = rho_ghost + DrhoDt*dt;
    } 
  


  	// singlephase or twophase PD-flash calculation
 	/*p_cal = p_ghost_update/1000.0; // kpa
  	rho_cal = rho_ghost_update/MW; // mol/L
  	PDFLSHdll(p_cal,rho_cal,x_cal,T_cal,rholiq_cal,rhovap_cal,xliq_cal,xvap_cal,q_cal,e_cal,h_cal,s_cal,cv_cal,cp_cal,w_cal,ierr,herr,255);
   
	T_ghost_update = T_cal;
  
  	if (q_cal >= 0.0 && q_cal <= 1.0) 
  	{
    		alpha1_ghost_update = q_cal*(rho_cal/rhovap_cal);
  	}
  	else if (q_cal > 1.0)
  	{
    		alpha1_ghost_update = 1.0;
  	}
  	else 
  	{
    		alpha1_ghost_update = 0.0;
  	}

  	Info << "p_ghost_update = " << p_ghost_update << nl << endl;	
	Info << "U_ghost_update = " << U_ghost_update << nl << endl;
	Info << "rho_ghost_update = " << rho_ghost_update << nl << endl;
	Info << "T_ghost_update = " << T_ghost_update << nl << endl;
  	Info << "alpha1_ghost_update = " << alpha1_ghost_update << nl << endl;*/
  
  //-----------------------------------------------------------
