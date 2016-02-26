/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    myTwoPhaseEulerFoam.C

Description
    currently a copy of twoPhaseEulerFoam only

    Solver for a system of 2 compressible fluid phases with one phase
    dispersed, e.g. gas bubbles in a liquid including heat-transfer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoPhaseSystem.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "IOMRFZoneList.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"

#include "phaseChangeTwoPhaseMixture.H"

#define REFPROP_IMPLEMENTATION
#include "REFPROP_lib.h"
#undef REFPROP_IMPLEMENTATION
//#include "myWaveTransmissiveFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createMRFZones.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "CourantNos.H"
    #include "setInitialDeltaT.H"
		
    //for testing only with a constant interpahse mass transfer term 	    	
    dimensionedScalar gamma_LV
 	("gamma_LV",
	 dimensionSet (1,-3,-1,0,0,0,0),
	 scalar(0.000));
	
    dimensionedScalar gamma_VL
 	("gamma_LV",
	 dimensionSet (1,-3,-1,0,0,0,0),
	 scalar(0.000));

    //constant and constant vectors
    vector unity(1,0,0);		
    
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
		celld = (double) (celli-(celln-100))/100;  
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



    //-------------------------------------------------------------------
    // Load the refprop thermodynamic libaray 
    //-------------------------------------------------------------------
    std::string err;
    bool loaded_REFPROP = load_REFPROP(err);
    long ierr = 0, nc = 1; 
    char herr[255], hfld[] = "CO2.FLD", hhmx[] = "HMX.BNC", href[] = "DEF";

    SETUPdll(nc,hfld,hhmx,href,ierr,herr,10000,255,3,255);

    Info<< "error massage:" << herr << nl << endl; 
    //-------------------------------------------------------------------	

    
    
    
    pimpleControl pimple(mesh);
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNos.H"
        #include "setDeltaT.H"

	//-----------------------------------------------------------
	// sonic/supersonic outflow boundary condition implementation
	// 
	// In 0/p file, sideRight boundary condition must be of 
	// fixed value type 
	//-----------------------------------------------------------		
	label patchID = mesh.boundaryMesh().findPatchID("sideRight");
	//p.boundaryField()[patchID] == scalar(1000000);
	const scalarField& p_test = p.boundaryField()[patchID];
	    
	// make a reference to cell adjacent to boundary patch 
	scalar p_interior = p[(celln)];
	scalar U_interior = mag(U[(celln)]);
  scalar rho_interior = alpha1[(celln)]*rho1[(celln)]+alpha2[(celln)]*rho2[(celln)];

	//Info << "pressure = " << p_interior << nl << endl;	
	//Info << "velocity = " << U_interior << nl << endl;
	//Info << "density = " << rho_interior << nl << endl;	
  
  
  // obtaining the ghost cell information 
	scalar U_ghost = mag(U.boundaryField()[patchID][0]);
  scalar p_ghost = p.boundaryField()[patchID][0];
	scalar rho_ghost = alpha1.boundaryField()[patchID][0]*rho1.boundaryField()[patchID][0]+alpha2.boundaryField()[patchID][0]*rho2.boundaryField()[patchID][0];
	scalar rho1_ghost = rho1.boundaryField()[patchID][0];
  scalar rho2_ghost = rho2.boundaryField()[patchID][0];
  scalar T1_ghost = thermo1.T().boundaryField()[patchID][0]; 
  scalar T2_ghost = thermo2.T().boundaryField()[patchID][0];
  scalar alpha1_ghost = alpha1.boundaryField()[patchID][0];
  scalar alpha2_ghost = alpha2.boundaryField()[patchID][0];

	
  Info << "pressure = " << p_ghost << nl << endl;	
	Info << "velocity = " << U_ghost << nl << endl;
	Info << "density = " << rho_ghost << nl << endl;	

  //
  // obtaining the speed of sound      
	//
	// Thermodynamic calculations
	
	long ierr = 0;
	char herr[255];
	double x_cal[] = {1.0}, xliq_cal[] = {0.0}, xvap_cal[] = {0.0}, q_cal = NULL,T_cal = NULL, rho_cal = NULL, rholiq_cal=NULL, rhovap_cal=NULL, p_cal=NULL, h_cal = NULL, s_cal = NULL, e_cal = NULL, cp_cal = NULL, cv_cal = NULL, w_cal = NULL;
	double hjt = NULL, dummy =NULL;	
  double MW = 44.01; // kg/kmol

  T_cal = T1_ghost, rho_cal = rho1_ghost/MW;
  THERM0dll(T_cal,rho_cal,x_cal,p_cal,e_cal,h_cal,s_cal,cv_cal,cp_cal,w_cal,hjt,dummy);
	   
  scalar sos1 = w_cal;

  T_cal = T2_ghost, rho_cal = rho2_ghost/MW;
  THERM0dll(T_cal,rho_cal,x_cal,p_cal,e_cal,h_cal,s_cal,cv_cal,cp_cal,w_cal,hjt,dummy);

	scalar sos2 = w_cal;
  // bulk speed of sound
  scalar sos = Foam::sqrt(1.0/(alpha1_ghost/(sos1*sos1)+alpha2_ghost/(sos2*sos2))); 
  //scalar sos = sqrt(static_cast<double>(sos_sqr));

  Info << "sound speed =" << sos << nl << endl;  
	
  
	//----------------------------------------------------------			
	
  
  // obtaining current timestep deltaT
  scalar dt = runTime.deltaT().value();
		
  // wave amplitude calculation
	if (runTime.time().value() == 0.0)
	{	
	    U_ghost = scalar(200);
	    p_ghost = scalar(1500000);
	    rho_ghost = scalar(20.19);
	}

  scalar p_ghost_update, U_ghost_update, rho_ghost_update;

  if (U_ghost > sos)
  {
	  scalar wave1 = 0;
	  scalar wave2 = (U_ghost/xDimDim[1])*(p_ghost - p_interior-sos*sos*(rho_ghost-rho_interior));
	  scalar wave3 = ((U_ghost+sos)/xDimDim[1])*(p_ghost-p_interior+rho_ghost*sos*(U_ghost-U_interior));
	  
    scalar DpDt = -0.5*(wave3+wave1);
	  scalar DuDt = -0.5/(rho_ghost*sos)*(wave3-wave1);
	  scalar DrhoDt = 1.0/(sos*sos)*(DpDt+wave2);

	  p_ghost_update = p_ghost + DpDt*dt;
	  U_ghost_update = U_ghost + DuDt*dt;
	  rho_ghost_update = rho_ghost + DrhoDt*dt;
  }
  else if (U_ghost <= sos)
  {
	  scalar wave2 = 0.0; //(U_ghost/xDimDim[1])*(p_ghost - p_interior-sos*sos*(rho_ghost-rho_interior));
	  scalar wave3 = ((U_ghost+sos)/xDimDim[1])*(p_ghost-p_interior+rho_ghost*sos*(U_ghost-U_interior));
    scalar wave1 = -wave3;
  
    
	  scalar DpDt = -0.5*(wave3+wave1);
	  scalar DuDt = -0.5/(rho_ghost*sos)*(wave3-wave1);
	  scalar DrhoDt = 1.0/(sos*sos)*(DpDt+wave2);

  	p_ghost_update = scalar(800000); //p_ghost + DpDt*dt;
	  U_ghost_update = U_ghost + DuDt*dt;
	  rho_ghost_update = rho_ghost + DrhoDt*dt;
  
  } 
  


  // singlephase or twophase PD-flash calculation
  p_cal = p_ghost_update/1000.0; // kpa
  rho_cal = rho_ghost_update/MW; // mol/L
  PDFLSHdll(p_cal,rho_cal,x_cal,T_cal,rholiq_cal,rhovap_cal,xliq_cal,xvap_cal,q_cal,e_cal,h_cal,s_cal,cv_cal,cp_cal,w_cal,ierr,herr,255);
   
  scalar T_ghost_update = T_cal;
  scalar alpha1_ghost_update;
  
  //if (q_cal >= 0.0 && q_cal <= 1.0) 
  //{
    //alpha1_ghost_update = q_cal*(rho_cal/rhovap_cal);
  //}
  //else if (q_cal > 1.0)
  //{
    alpha1_ghost_update = 1.0;
  //}
  //else 
  //{
    //alpha1_ghost_update = 0.0;
  //}


  Info << "p_ghost_update = " << p_ghost_update << nl << endl;	
	Info << "U_ghost_update = " << U_ghost_update << nl << endl;
	Info << "rho_ghost_update = " << rho_ghost_update << nl << endl;
	Info << "T_ghost_update = " << T_ghost_update << nl << endl;
  Info << "alpha1_ghost_update = " << alpha1_ghost_update << nl << endl;
  
  //-----------------------------------------------------------
        
	// runtime time output    
	runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fluid.solve();
            fluid.correct();

            // update boundary conditions
	          p.boundaryField()[patchID] == p_ghost_update;
            U1.boundaryField()[patchID] == vector(U_ghost_update,0,0);
            U2.boundaryField()[patchID] == vector(0,0,0); //vector(U_ghost_update,0,0);
            thermo1.T().boundaryField()[patchID] == T_ghost_update;
            thermo2.T().boundaryField()[patchID] == T_ghost_update;

	          U_bulk = mag(alpha1*U1+alpha2*U2);
	          rho_bulk = alpha1*rho1+alpha2*rho2;						
            psi_bulk =1.0/(alpha1/thermo1.psi()+alpha2/thermo2.psi());		 
	    
            volScalarField contErr1
            (
                fvc::ddt(alpha1, rho1) + fvc::div(alphaRhoPhi1)
              - (fvOptions(alpha1, rho1)&rho1)                       //+gamma_LV-gamma_VL // fvOptions are the runtime semi-implicit source term 
              + alpha1*rho1*mag(U1)*areaSource	
	          );

            volScalarField contErr2
            (
                fvc::ddt(alpha2, rho2) + fvc::div(alphaRhoPhi2)
               - (fvOptions(alpha2, rho2)&rho2)                    //-gamma_LV+gamma_VL // 
               + alpha2*rho2*mag(U2)*areaSource	 
	          );

			
            #include "UEqns.H"

	          U_bulk = mag(alpha1*U1+alpha2*U2);                     // update velocity field                   				
       
            #include "EEqns.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            #include "DDtU.H"

            if (pimple.turbCorr())
            {
                fluid.correctTurbulence();
            }
        }

        #include "write.H"
	
	//const volScalarField& test = alpha1_.db().lookupObject<volScalarField>("flowAreaGrad");


        Info<< "flowArea=" << flowArea[950] << nl << endl; 
   	//Info<< "flowAreaGrad=" << flowAreaGrad[950] << nl << endl; 
        //Info<< "areaSource=" << areaSource[950] << nl << endl;


// loop over all cells:
	//forAll(mesh.C(), cellI) 
	//{
		//Info << "******* CellID: " << cellI << "*******"<< endl;

		//Getting list of all faces of current cell
		//const labelList& faces = mesh.cells()[cellI];

		//loop over all faces of current cell
		//forAll( faces, faceI )
		//{   
			//if (mesh.isInternalFace(faces[faceI]))
			//{   
				//Info << "internal faceI: " << faceI << "    mesh.Sf()[faceI]: " << -mesh.Sf()[faceI] << "    mesh.magSf()[faceI]: " << mesh.magSf()[faceI] << endl;
			//}   
			//else
			//{   
				//Info << "boundary faceI: " << faceI << "    mesh.Sf()[faceI]: " << -mesh.Sf()[faceI] << "    mesh.magSf()[faceI]: " << mesh.magSf()[faceI] << endl;
			//}   

		//} //move on to next face  
		//Info << " " << endl;
	//}//move on to next cell 


	Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
