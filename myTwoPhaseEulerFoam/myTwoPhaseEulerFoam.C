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

    // pipeline dimension 
    dimensionedScalar Dh
 	("Dh",
	 dimensionSet (0,1,0,0,0,0,0),
	 scalar(0.233));
   
    // friction factor
 	  scalar frictionFactor = 0.005;

    //my volScalarField declaration 
    #include "myVolScalar.H"

    //quasi-one-dimensional flow setup
    #include "changeArea.H"

    //load refprop thermodynamic library 			
    #include "refpropLibLoading.H"	

    //dummy vector  
    vector unity(1,0,0);

    //start of the loop       
    pimpleControl pimple(mesh);
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNos.H"
        #include "setDeltaT.H"

	// interface mass and heat transfer 
	#include "massAndEnergyTransfer.H"
	
	// update boundary conditions (NSCBC)
	#include "NSCBC.H"
	
	// runtime time output    
	runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        { 

            
            fluid.solve();
            fluid.correct();

	    // interface mass and heat transfer 
            //#include "massAndEnergyTransfer.H"

            
	    // update boundary conditions
	    p.boundaryField()[patchID] == p_ghost_update;
            U1.boundaryField()[patchID] == vector(U_ghost_update,0,0);
            U2.boundaryField()[patchID] == vector(U_ghost_update,0,0);
            //thermo1.T().boundaryField()[patchID] == T_ghost_update;
            //thermo2.T().boundaryField()[patchID] == T_ghost_update;           
	    //alpha1.boundaryField()[patchID] == alpha1_ghost_update;
	    //alpha2.boundaryField()[patchID] == 1.0 - alpha1_ghost_update;


            U_bulk = mag(alpha1*U1+alpha2*U2);
	    rho_bulk = alpha1*rho1+alpha2*rho2;						
            psi_bulk =1.0/(alpha1/thermo1.psi()+alpha2/thermo2.psi());		 
	    
            volScalarField contErr1
            (
                fvc::ddt(alpha1, rho1) + fvc::div(alphaRhoPhi1)
              - (fvOptions(alpha1, rho1)&rho1)                      // fvOptions are the runtime semi-implicit source term 
              + alpha1*rho1*mag(U1)*areaSource
	      - gammaV
	          );

            volScalarField contErr2
            (
                fvc::ddt(alpha2, rho2) + fvc::div(alphaRhoPhi2)
               - (fvOptions(alpha2, rho2)&rho2)                     
               + alpha2*rho2*mag(U2)*areaSource
	       - gammaL
	          );
            
            // update friction source term  
	    Fv = - scalar(2)*frictionFactor*alpha1*rho1*mag(U1)*mag(U1)/Dh;  
            Fl = - scalar(2)*frictionFactor*alpha2*rho2*mag(U2)*mag(U2)/Dh;
            
            #include "UEqns.H"

            // update friction source term for energy balance 
	    U_bulk = mag(alpha1*U1+alpha2*U2);                    				
            FvU_bulk = U_bulk*Fv;
            FlU_bulk = U_bulk*Fl;

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

	//loop over all cells:
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
