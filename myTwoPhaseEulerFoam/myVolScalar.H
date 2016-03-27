    // declaration of cell-to-cell length
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
    // declaratino of the cell centre variable 
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
   

    // declaration of the gradient of area along flow direction
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
    
    // declaration of the area part of the source term
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
    
    // declaration of  bulk density 
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


    // declaration of transporting variables 
    volScalarField qV
    (
     IOobject
     (
 	"qV",
 	runTime.timeName(),
 	mesh
     ),
	mesh,
 	dimensionedScalar("0",dimensionSet (1,-1,-3,0,0),0)
    );


    volScalarField qL
    (
     IOobject
     (
 	"qL",
 	runTime.timeName(),
 	mesh
     ),
	mesh,
 	dimensionedScalar("0",dimensionSet (1,-1,-3,0,0),0)
    );


    volScalarField gammaV
    (
     IOobject
     (
 	"gammaV",
 	runTime.timeName(),
 	mesh
     ),
	mesh,
 	dimensionedScalar("0",dimensionSet (1,-3,-1,0,0),0)
    );


   volScalarField gammaL
   (
     IOobject
     (
 	"gammaL",
 	runTime.timeName(),
 	mesh
     ),
	mesh,
 	dimensionedScalar("0",dimensionSet (1,-3,-1,0,0),0)
   );


   volScalarField gammaVu
   (
     IOobject
     (
 	"gammaVu",
 	runTime.timeName(),
 	mesh
     ),
	mesh,
 	dimensionedScalar("0",dimensionSet (1,-2,-2,0,0),0)
   );


   volScalarField gammaLu
   (
     IOobject
     (
 	"gammaLu",
 	runTime.timeName(),
 	mesh
     ),
	mesh,
 	dimensionedScalar("0",dimensionSet (1,-2,-2,0,0),0)
   );

   // declaration of the friction souce terms 	
   volScalarField Fv
   (
     IOobject
     (
 	"Fv",
 	runTime.timeName(),
 	mesh
     ),
	mesh,
 	dimensionedScalar("0",dimensionSet (1,-2,-2,0,0),0)
   );

   volScalarField Fl
   (
     IOobject
     (
 	"Fl",
 	runTime.timeName(),
 	mesh
     ),
	mesh,
 	dimensionedScalar("0",dimensionSet (1,-2,-2,0,0),0)
   );


   volScalarField FvU_bulk
   (
     IOobject
     (
 	"FvU_bulk",
 	runTime.timeName(),
 	mesh
     ),
	mesh,
 	dimensionedScalar("0",dimensionSet (1,-1,-3,0,0),0)
   );


   volScalarField FlU_bulk
   (
     IOobject
     (
 	"FlU_bulk",
 	runTime.timeName(),
 	mesh
     ),
	mesh,
 	dimensionedScalar("0",dimensionSet (1,-1,-3,0,0),0)
   );
