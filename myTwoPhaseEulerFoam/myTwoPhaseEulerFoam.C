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
    //label patchID = mesh.boundaryMesh().findPatchID("sideLeft");//locate particular patch ID

    //Info<< "patchID" << patchID << nl << endl; 
    
    //const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
    //patch().magSf()[patchID] = patch().magSf()[patchID]*scalar(0.1);				

    // Load the shared library
    //-------------------------------------------------------------------
    //std::string err;
    //bool loaded_REFPROP = load_REFPROP(err);
    //long ierr = 0, nc = 1; 
    //char herr[255], hfld[] = "CO2.FLD", hhmx[] = "HMX.BNC", href[] = "DEF";

    //SETUPdll(nc,hfld,hhmx,href,ierr,herr,10000,255,3,255);

    //Info<< "error massage:" << herr << nl << endl; 
    //-------------------------------------------------------------------	

    pimpleControl pimple(mesh);
      
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNos.H"
        #include "setDeltaT.H"


	//Thermodynamic calculations
	//{
		//long ierr = 0;
		//char herr[255];
		//double z[] = {1.0}, x[] = {1.0}, y[] = {1.0}, T= 300, p = 1010.325, d = NULL, dl = NULL, dv = NULL, h = NULL, s = NULL, u = NULL, cp = NULL, cv = NULL, q = NULL, w = NULL;
		//TPFLSHdll(T, p, z, d, dl, dv, x, y, h,s,u,cp,cv,w,q,ierr,herr,255);
	
	
       	//Info<< "enthalpy = " << w << nl << endl; 
	
	//}


        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fluid.solve();
            fluid.correct();

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
