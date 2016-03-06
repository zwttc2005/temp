//-----------------------------------------------------------------------------\\
// Interphase energy transfer calculation based on Barcia-Cascales et al., 2007
//
// It is a linear relaxation equation (relaxed from saturation)
//
//-----------------------------------------------------------------------------\\

// define problem specific constant
scalar Ctou = 1.0, tou = 0.00001;

// hSatV, hSatL
// hSatV = -0.000508*Foam::pow (T,4)+0.402133*Foam::pow(T,3)-121.645905*Foam::pow(T,2)+17114.230257*T-537903.887493;
//
// hSatL = 0.000491*Foam::pow(T,4)-0.397733*Foam::pow(T,3)+120.103770*Foam::pow(T,2)-14045.518453*T+449986.580984;

volScalarField& Tvap = thermo1.T();
volScalarField& Tliq = thermo2.T();
volScalarField& hVap = thermo1.he();
volScalarField& hLiq = thermo2.he();

// local variables
scalar TvapS, TliqS, pS, hVapS, hLiqS, rhoVapS, rhoLiqS, qVs, qLs, gammaLs, gammaVs;
scalar pSatV, pSatL, hSatV, hSatL;

int i;

for(i = 0; i < celln+1; i++)
{
	TvapS = Tvap[i];//.value();
	TliqS = Tliq[i];//.value();
	
	// equilibrium relation
	
	hSatV = -0.000508*Foam::pow (TvapS,4)+0.402133*Foam::pow(TvapS,3)-121.645905*Foam::pow(TvapS,2)+17114.230257*TvapS-537903.887493;
	hSatL = 0.000491*Foam::pow(TliqS,4)-0.397733*Foam::pow(TliqS,3)+120.103770*Foam::pow(TliqS,2)-14045.518453*TliqS+449986.580984;
	
	pSatV = 2.600717*Foam::pow(TvapS,3)-1288.694897*Foam::pow(TvapS,2)+215032.266882*TvapS-12071874.471421;
	pSatL = 2.600717*Foam::pow(TliqS,3)-1288.694897*Foam::pow(TliqS,2)+215032.266882*TliqS-12071874.471421;

	// inheriet dimensioned scalar to scalar
	hVapS = hVap[i];
	hLiqS = hLiq[i];
	rhoVapS = rho1[i];
	rhoLiqS = rho2[i];
	pS = p[i];


	if (alpha1[i] > 0.0 && alpha2[i] > 0.0) 
	{	
		qV[i] = (Ctou/tou)*alpha1[i]*alpha2[i]*(hSatV-hVapS);
		qL[i] = (Ctou/tou)*alpha1[i]*alpha2[i]*(hSatL-hLiqS);

		// mass transfer calculation
		qVs = qV[i];
		qLs = qL[i];
	
		gammaV[i] = -(qVs+qLs)/(hSatV-hSatL);
		gammaL[i] = -gammaV[i];
	
		if (gammaVs >= 0.0)
		{
			//qV[i] = gammaVs*hSatV;
			//qL[i] = gammaLs*hSatV;
		
			gammaVu[i] = gammaV[i]*mag(U1[i]);
			gammaLu[i] = gammaL[i]*mag(U1[i]);
	
		}
		else if (gammaVs < 0.0)
		{
			//qV[i] = gammaVs*hSatL;
			//qL[i] = gammaLs*hSatL;

			gammaVu[i] = gammaV[i]*mag(U2[i]);
			gammaLu[i] = gammaL[i]*mag(U2[i]);
		}
	}
	else
	{ 
		// private doubles
		scalar mDotV, mDotL;

		// mass trasnfer

		if (i == 999) Info << "T_999 " << TliqS << nl << endl;
		if (i == 999) Info << "Delta " << pSatL << nl << endl;
		if (i == 999) Info << "liquid " << alpha2[i] << nl << endl;

		if (TvapS < 300.0 && alpha1[i] > 0.0 && pS > pSatV && hVapS < hSatV)
		{
			mDotV = (Ctou/0.001)*(alpha1[i]*rhoVapS*((pS-pSatV)/pSatV)); // condensation	
			
			Info << "massflux" << mDotV << nl << endl;
		}  
		else
		{
			mDotV = 0.0;
		}

		if (TliqS < 300.0 && alpha2[i] > 0.0 && pSatL > pS && hLiqS > hSatL) 
		{
			mDotL = (Ctou/0.001)*(alpha2[i]*rhoLiqS*((pSatL-pS)/pSatL)); // evaporation 
		}
		else
		{
			mDotL = 0.0;
		}
		
		gammaV[i] = mDotL-mDotV;
		gammaL[i] = - gammaV[i];
		gammaVs = gammaV[i];
		gammaLs = gammaL[i]; 

		if (gammaVs >= 0.0)
		{
			qV[i] = gammaVs*hSatV;
			qL[i] = gammaLs*hSatV;
		
			gammaVu[i] = gammaV[i]*mag(U1[i]);
			gammaLu[i] = gammaL[i]*mag(U1[i]);
	
		}
		else if (gammaVs < 0.0)
		{
			qV[i] = gammaVs*hSatL;
			qL[i] = gammaLs*hSatL;

			gammaVu[i] = gammaV[i]*mag(U2[i]);
			gammaLu[i] = gammaL[i]*mag(U2[i]);
		}

	}		

	
		//qL[i] = scalar(0);
		//gammaL[i] = scalar(0);
		//gammaLu[i] = scalar(0);
}
    


//{	
	//double Tgas,pgas,rhogas,e=NULL,h=NULL,s=NULL,Cv=NULL,Cp=NULL,w=NULL,hjt=NULL,dummy=NULL;   	
	//double xgas[]={1.0};
    	//long kph = -2, kguess = 0; // integer required for TP-calculation
	
	//Tgas = Tvap[999];
    	//double pSatV = 2.600717*Foam::pow(Tgas,3)-1288.694897*Foam::pow(Tgas,2)+215032.266882*Tgas-12071874.471421;
    	
	//pgas = pSatV/1000.0;

	//TPRHOdll(Tgas,pgas,xgas,kph,kguess,rhogas,ierr,herr,255);	

	//Tgas = 260.0, rhogas = 2.0; 
		//Info << "psat = " << pSatV << nl << endl;
		//Info << "pgas_sat = " << pgas << nl << endl;		
		//Info << "rhosat = " << rhogas*44.01 << nl << endl;
		//Info << "Tgas = " << Tgas << nl << endl;
		
	//THERMdll(Tgas,rhogas,xgas,pgas,e,h,s,Cv,Cp,w,hjt);
    

    	//if (p < 2500000)
    	//{	 		
		//Info << "psat = " << pSatV << nl << endl;
		//Info << "pgas_sat = " << pgas << nl << endl;		
		//Info << "rhosat = " << rhogas*44.01 << nl << endl;
		//Info << "Tgas = " << Tgas << nl << endl;
		//Info << "hsat = " << (h*1000.0/44.01) << nl << endl;
		//Info << "hVapS = " << hVap[999] << nl << endl;
 
    		//Info << " Tgas = " << Tgas << nl << endl;
		//Info << "rhosat = " << rhogas << nl << endl;
		//Info << "hsat = " << (h*1000.0/44.01) << nl << endl; 
    	//}
//}

Info << "gammaL = " << gammaL[999] << nl << endl;
Info << "qL = " << qL[999] << nl << endl;
Info << "hLiqS = " << hLiq[999] << nl << endl;
Info << "hSatL = " << hSatL << nl << endl;


