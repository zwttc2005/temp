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


    // define transporting variables 
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
