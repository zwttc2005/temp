//-----------------------------------------------------------------------------\\
// Interphase energy transfer calculation based on Barcia-Cascales et al., 2007
//
// It is a linear relaxation equation (relaxed from saturation)
//
//-----------------------------------------------------------------------------\\

// define problem specific constants
   scalar Ctou = 1.0, tou = 0.00005; //0.00001;

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
scalar pSatV, pSatL, hSatV, hSatL,hSatVatTliqS,hSatLatTvapS;

int i;

for(i = 0; i < celln+1; i++)
{
	TvapS = Tvap[i];//.value();
	TliqS = Tliq[i];//.value();
	
	// equilibrium relation
	
	//hSatV = -0.000499*Foam::pow (TvapS,4)+0.391180*Foam::pow(TvapS,3)-116.809829*Foam::pow(TvapS,2)+16284.795610*TvapS-493979.751650;
	//hSatL = 0.000761*Foam::pow(TliqS,4)-0.622377*Foam::pow(TliqS,3)+191.786710*Foam::pow(TliqS,2)-24634.683426*TliqS+1069742.444044;
	
	//pSatV = 2.600717*Foam::pow(TvapS,3)-1288.694897*Foam::pow(TvapS,2)+215032.266882*TvapS-12071874.471421;
	//pSatL = 2.600717*Foam::pow(TliqS,3)-1288.694897*Foam::pow(TliqS,2)+215032.266882*TliqS-12071874.471421;


	hSatV = (-0.00033568697770647300*Foam::pow (TvapS,4)+0.22286341085225600000*Foam::pow(TvapS,3)-56.34175837976970000000*Foam::pow(TvapS,2)+7091.43647987876000000000*TvapS+11004.31020779390000000000)*
		       (1.0-(0.5+0.5*Foam::tanh((TvapS-216.6)/5.0)))+
		       (-0.00000274435901514991*Foam::pow (TvapS,6)+0.00416370607053285000*Foam::pow(TvapS,5)-2.62812735352465000000*Foam::pow(TvapS,4)
		       +883.27046130716200000000*Foam::pow(TvapS,3)-166692.79608873600000000000*Foam::pow(TvapS,2)+16749035.85359570000000000000*TvapS-699609750.01779100000000000000)*
		       (0.5+0.5*Foam::tanh((TvapS-216.6)/5.0));

	hSatL = (0.00005160214658772590*Foam::pow (TliqS,4)-0.02405779531282800000*Foam::pow(TliqS,3)+4.64613689715136000000*Foam::pow(TliqS,2)+1060.54667046420000000000*TliqS-236588.78011277400000000000)*
		(1.0-(0.5+0.5*Foam::tanh((TliqS-216.6)/5.0)))+
		(0.00472233079677409000 *Foam::pow (TliqS,4)-4.66454799246256000000*Foam::pow(TliqS,3)+1731.11661128817000000000*Foam::pow(TliqS,2)-283951.38280328500000000000*TliqS+17375579.55332360000000000000)*
		(0.5+0.5*Foam::tanh((TliqS-216.6)/5.0));
	

	pSatV = (2.00185689022271000000*Foam::pow (TvapS,3)-946.33910600661800000000*Foam::pow(TvapS,2)+150059.16194839700000000000*TvapS-7978874.90064541000000000000)*
		(1.0-(0.5+0.5*Foam::tanh((TvapS-216.6)/5.0)))+
		(3.66785811262889000000*Foam::pow (TvapS,3)-2075.51035670498000000000*Foam::pow(TvapS,2)+406181.92100634200000000000*TvapS-27366391.52576530000000000000)*
		(0.5+0.5*Foam::tanh((TvapS-216.6)/5.0));
	

	pSatL = (2.00185689022271000000*Foam::pow (TliqS,3)-946.33910600661800000000*Foam::pow(TliqS,2)+150059.16194839700000000000*TliqS-7978874.90064541000000000000)*
		(1.0-(0.5+0.5*Foam::tanh((TliqS-216.6)/5.0)))+
		(3.66785811262889000000*Foam::pow (TliqS,3)-2075.51035670498000000000*Foam::pow(TliqS,2)+406181.92100634200000000000*TliqS- 27366391.52576530000000000000)*
		(0.5+0.5*Foam::tanh((TliqS-216.6)/5.0));



	// specific enthalpy that is transported 
	hSatVatTliqS = (-0.00033568697770647300*Foam::pow (TliqS,4)+0.22286341085225600000*Foam::pow(TliqS,3)-56.34175837976970000000*Foam::pow(TliqS,2)+7091.43647987876000000000*TliqS+11004.31020779390000000000)*
		       (1.0-(0.5+0.5*Foam::tanh((TliqS-216.6)/5.0)))+
		       (-0.00000274435901514991*Foam::pow (TliqS,6)+0.00416370607053285000*Foam::pow(TliqS,5)-2.62812735352465000000*Foam::pow(TliqS,4)
		       +883.27046130716200000000*Foam::pow(TliqS,3)-166692.79608873600000000000*Foam::pow(TliqS,2)+16749035.85359570000000000000*TliqS-699609750.01779100000000000000)*
		       (0.5+0.5*Foam::tanh((TliqS-216.6)/5.0));




	hSatLatTvapS = (0.00005160214658772590*Foam::pow (TvapS,4)-0.02405779531282800000*Foam::pow(TvapS,3)+4.64613689715136000000*Foam::pow(TvapS,2)+1060.54667046420000000000*TvapS-236588.78011277400000000000)*
		(1.0-(0.5+0.5*Foam::tanh((TvapS-216.6)/5.0)))+
		(0.00472233079677409000 *Foam::pow (TvapS,4)-4.66454799246256000000*Foam::pow(TvapS,3)+1731.11661128817000000000*Foam::pow(TvapS,2)-283951.38280328500000000000*TvapS+17375579.55332360000000000000)*
		(0.5+0.5*Foam::tanh((TvapS-216.6)/5.0));





	// inheriet dimensioned scalar to scalar
	hVapS = hVap[i];
	hLiqS = hLiq[i];
	rhoVapS = rho1[i];
	rhoLiqS = rho2[i];
	pS = p[i];


	//if (alpha1[i] > 0.0 && alpha2[i] > 0.0) 
	{	
		if (hSatV > hVapS && pS > pSatV)
    		{
        		qV[i] = (Ctou/tou)*alpha1[i]*alpha2[i]*(hSatV-hVapS);
    		}
    		else 
    		{	
        		qV[i] = scalar(0);
    		}
    
    		if (hSatL < hLiqS && pS < pSatL)    
    		{    
        		qL[i] = (Ctou/tou)*alpha1[i]*alpha2[i]*(hSatL-hLiqS);
    		}
    		else
    		{
        		qL[i] = scalar(0);
    		}
    
    
    		// mass transfer calculation
		qVs = qV[i];
		qLs = qL[i];
	
		gammaV[i] = -(qVs+qLs)/(hSatV-hSatL);
		gammaL[i] = -gammaV[i];
	
		if (gammaVs >= 0.0)
		{
			gammaVu[i] = gammaV[i]*mag(U1[i]);
			gammaLu[i] = gammaL[i]*mag(U1[i]);
	
		}
		else if (gammaVs < 0.0)
		{

			gammaVu[i] = gammaV[i]*mag(U2[i]);
			gammaLu[i] = gammaL[i]*mag(U2[i]);
		}
	}
  //else
	
	/*{ 
		// private doubles
		scalar mDotV, mDotL;

		// mass trasnfer

		if (TvapS < 300.0 && alpha1[i] > 0.0 && pS > pSatV)                // hVapS < hSatV
		{
			mDotV = (Ctou/tou)*(alpha1[i]*rhoVapS*((pS-pSatV)/pSatV)); // condensation	
			
			//Info << "massflux" << mDotV << nl << endl;
		}  
		else
		{
			mDotV = 0.0;
		}

		if (TliqS < 300.0 && alpha2[i] > 0.0 && pSatL > pS)                // hLiqS > hSatL
		{
			mDotL = (Ctou/tou)*(alpha2[i]*rhoLiqS*((pSatL-pS)/pSatL)); // evaporation 
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
			qV[i] = gammaVs*hSatVatTliqS;
			qL[i] = gammaLs*hSatVatTliqS;
		
			gammaVu[i] = gammaV[i]*mag(U1[i]);
			gammaLu[i] = gammaL[i]*mag(U1[i]);
	
		}
		else if (gammaVs < 0.0)
		{
			qV[i] = gammaVs*hSatLatTvapS;
			qL[i] = gammaLs*hSatLatTvapS;

			gammaVu[i] = gammaV[i]*mag(U2[i]);
			gammaLu[i] = gammaL[i]*mag(U2[i]);
		}

	}*/		

	
	/*{ 
		// private doubles
		scalar mDotV, mDotL;

		// mass trasnfer

		//if (i == 999) Info << "T_999 " << TliqS << nl << endl;
		//if (i == 999) Info << "Delta " << pSatL << nl << endl;
		//if (i == 999) Info << "liquid " << alpha2[i] << nl << endl;

		if (TvapS < 300.0 && alpha1[i] > 0.0 && hVapS < hSatV)
		{
			mDotV = (Ctou/tou)*rho1[i]*alpha1[i]*(1.0-(hVapS-hSatLatTvapS)/(hSatV-hSatLatTvapS)); // condensation	
			
			//Info << "massflux" << mDotV << nl << endl;
		}  
		else
		{
			mDotV = 0.0;
		}

		if (TliqS < 300.0 && alpha2[i] > 0.0 && hLiqS > hSatL)
		{
			// saturated vapour phase enthalpy at liquid phase temperature 
	    		//hSatVatTliqS = -0.000499*Foam::pow (TliqS,4)+0.391180*Foam::pow(TliqS,3)-116.809829*Foam::pow(TliqS,2)+16284.795610*TliqS-493979.751650;
      		
			mDotL = (Ctou/tou)*rho2[i]*alpha2[i]*((hLiqS-hSatL)/(hSatVatTliqS-hSatL)); // evaporation 
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
			qV[i] = gammaVs*hSatVatTliqS;
			
			//Info << "hSatVatTliqS" << hSatVatTliqS<<nl<<endl;
			//Info << "hSatLatTvapS" << hSatLatTvapS<<nl<<endl;
			//Info << "hSatV" << hSatV<<nl<<endl;	
	
			//Info << "hSatL" << hSatL<<nl<<endl;
			//Info << "pSatV" << pSatV<<nl<<endl;
			//Info << "pSatL" << pSatL<<nl<<endl;

			qL[i] = gammaLs*hSatVatTliqS;
		
			gammaVu[i] = gammaV[i]*mag(U1[i]);
			gammaLu[i] = gammaL[i]*mag(U1[i]);
	
		}
		else if (gammaVs < 0.0)
		{
			qV[i] = gammaVs*hSatLatTvapS;
			qL[i] = gammaLs*hSatLatTvapS;

			gammaVu[i] = gammaV[i]*mag(U2[i]);
			gammaLu[i] = gammaL[i]*mag(U2[i]);
		}

	}*/


	/*qL[i] = scalar(0);
	gammaL[i] = scalar(0);
	gammaLu[i] = scalar(0);
	qV[i] = scalar(0);
	gammaV[i] = scalar(0);
	gammaVu[i] = scalar(0);*/		


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

Info << "gammaL = " << gammaL[(celln-10)] << nl << endl;
Info << "qL = " << qL[(celln-10)] << nl << endl;
Info << "hLiqS = " << hLiq[(celln-10)] << nl << endl;
//Info << "hSatL = " << hSatL << nl << endl;
