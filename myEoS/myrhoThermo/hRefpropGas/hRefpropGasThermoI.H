/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

// add refprop header here 
#define REFPROP_IMPLEMENTATION
#include "REFPROP_lib.h"
#undef REFPROP_IMPLEMENTATION

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class equationOfState>
inline Foam::hRefpropGasThermo<equationOfState>::hRefpropGasThermo
(
    const equationOfState& st,
    const scalar cp,
    const scalar hf
)
:
    equationOfState(st),
    Cp_(cp),
    Hf_(hf)
{}





// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::hRefpropGasThermo<equationOfState>::hRefpropGasThermo
(
    const word& name,
    const hRefpropGasThermo& ct
)
:
    equationOfState(name, ct),
    Cp_(ct.Cp_),
    Hf_(ct.Hf_)
{}


template<class equationOfState>
inline Foam::autoPtr<Foam::hRefpropGasThermo<equationOfState> >
Foam::hRefpropGasThermo<equationOfState>::clone() const
{
    return autoPtr<hRefpropGasThermo<equationOfState> >
    (
        new hRefpropGasThermo<equationOfState>(*this)
    );
}


template<class equationOfState>
inline Foam::autoPtr<Foam::hRefpropGasThermo<equationOfState> >
Foam::hRefpropGasThermo<equationOfState>::New(Istream& is)
{
    return autoPtr<hRefpropGasThermo<equationOfState> >
    (
        new hRefpropGasThermo<equationOfState>(is)
    );
}


template<class equationOfState>
inline Foam::autoPtr<Foam::hRefpropGasThermo<equationOfState> >
Foam::hRefpropGasThermo<equationOfState>::New(const dictionary& dict)
{
    return autoPtr<hRefpropGasThermo<equationOfState> >
    (
        new hRefpropGasThermo<equationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class EquationOfState>
inline Foam::scalar Foam::hRefpropGasThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class equationOfState>
inline Foam::scalar Foam::hRefpropGasThermo<equationOfState>::cp
(
    const scalar p,
    const scalar T

) const
{    
    scalar cp_;	 
    scalar rhogas;
        
    /*scalar Tc_ = 304.13;
    scalar Vc_ = 0.0021386/this->W();
    scalar Pc_ = 7377300;
    scalar Zc_ = Pc_*Vc_/(this->R()*Tc_);	
    scalar omega_ = 0.22394;*/


    // Load the shared library
    //-------------------------------------------------------------------
    
    static int flag = 0;
    long ierr = 0;
    char herr[255];

    if (flag == 0)
    {
    	std::string err;
    	bool loaded_REFPROP = load_REFPROP(err);
    	long nc = 1, n_PREOS = 2; 
    	char hfld[] = "CO2.FLD", hhmx[] = "HMX.BNC", href[] = "DEF";

    	SETUPdll(nc,hfld,hhmx,href,ierr,herr,10000,255,3,255);
        
	// use PR-EoS
	PREOSdll(n_PREOS);
    	if (ierr != 0)
    	{
    		Info<< "error massage:" << herr << nl << endl;  
   	}
    	
	flag++;
    }
    //-------------------------------------------------------------------		
    /*{
    	scalar a = 0.45724*sqr(this->R())*sqr(Tc_)/Pc_;
    	scalar b = 0.07780*this->R()*Tc_/Pc_;
    	scalar Tr = T/Tc_;

    	scalar alpha =
        	sqr
        	(
          	1.0
          	+ (0.37464 + 1.54226*omega_- 0.26992*sqr(omega_))
          	* (1.0 - sqrt(Tr))
        	);

    	scalar B = b*p/(this->R()*T);
    	scalar A = a*alpha*p/sqr(this->R()*T);

    	scalar a2 = B - 1.0;
    	scalar a1 = A - 2.0*B - 3.0*sqr(B);
    	scalar a0 = -A*B + sqr(B) + pow3(B);

    	scalar Q = (3.0*a1 - a2*a2)/9.0;
    	scalar Rl = (9.0*a2*a1 - 27.0*a0 - 2.0*a2*a2*a2)/54;

    	scalar Q3 = Q*Q*Q;
    	scalar D = Q3 + Rl*Rl;

    	scalar root = -1.0;
	
    	if (D <= 0.0)
    	{
        	scalar th = ::acos(Rl/sqrt(-Q3));
        	scalar qm = 2.0*sqrt(-Q);
       	 	scalar r1 = qm*cos(th/3.0) - a2/3.0;
        	scalar r2 = qm*cos((th + 2.0*constant::mathematical::pi)/3.0) - a2/3.0;
        	scalar r3 = qm*cos((th + 4.0*constant::mathematical::pi)/3.0) - a2/3.0;

        	root = max(r1, max(r2, r3));
		rhogas = p/(root*this->R()*T);
    	}
    	else
   	{
        	//one root is real
        	scalar D05 = sqrt(D);
        	scalar S = pow(Rl + D05, 1.0/3.0);
        	scalar Tl = 0;
        	if (D05 > Rl)
       	 	{
            		Tl = -pow(mag(Rl - D05), 1.0/3.0);
        	}
        	else
        	{
            		Tl = pow(Rl - D05, 1.0/3.0);
        	}
        
		root = S + Tl - a2/3.0;
		rhogas = p/(root*this->R()*T);
    	}
    }*/	
    //-------------------------------------------------------------------------
    
    {
    	double Tgas,pgas,e=NULL,h=NULL,s=NULL,Cv=NULL,Cp=NULL,w=NULL,hjt=NULL;
    	double xgas[]={1.0};
    	//long kph = -2, kguess = 0; // integer required for TP-calculation

    	Tgas = T;        //K
	rhogas = this->rho(p,T)/this->W();
    	//pgas = p/1000.0; //kpa

        
	//TPRHOdll(Tgas,pgas,xgas,kph,kguess,rhogas,ierr,herr,255);	
	
    	//THERMdll(Tgas,rhogas,xgas,pgas,e,h,s,Cv,Cp,w,hjt);
	
	CVCPdll(Tgas,rhogas,xgas,Cv,Cp);

    	if (p < 2500000)
    	{
      		//Info << "internal density = " << rhogas << nl << endl;
    		//Info << "internal temperature = " << Tgas << nl << endl;
      		//Info << "hs = " << (h*1000.0/44.01) << nl << endl;  
    		//Info << "Cp*T = " << Cp_*T << nl << endl;
    	}
	
	cp_ = Cp*1000.0;
    }
    
    return (cp_);
}


template<class equationOfState>
inline Foam::scalar Foam::hRefpropGasThermo<equationOfState>::ha
(
    const scalar p, const scalar T
) const
{
    return Cp_*T + Hf_;
}


template<class equationOfState>
inline Foam::scalar Foam::hRefpropGasThermo<equationOfState>::hs
(
    const scalar p, const scalar T
) const
{
    scalar hs_;
    scalar rhogas;	

    //scalar Tc_ = 304.13;
    //scalar Vc_ = 0.0021386/this->W();
    //scalar Pc_ = 7377300;
    //scalar Zc_ = Pc_*Vc_/(this->R()*Tc_);	
    //scalar omega_ = 0.22394;

    // Load the shared library
    //-------------------------------------------------------------------
    
    static int flag = 0;
    long ierr = 0;
    char herr[255];

    if (flag == 0)
    {
    	std::string err;
    	bool loaded_REFPROP = load_REFPROP(err);
    	long nc = 1, n_PREOS = 2; 
    	char hfld[] = "CO2.FLD", hhmx[] = "HMX.BNC", href[] = "DEF";

    	SETUPdll(nc,hfld,hhmx,href,ierr,herr,10000,255,3,255);
        
	// use PR-EoS
	PREOSdll(n_PREOS);
    	if (ierr != 0)
    	{
    		Info<< "error massage:" << herr << nl << endl;  
   	}
    	
	flag++;
    }
    //-------------------------------------------------------------------		
    /*{
    	scalar a = 0.45724*sqr(this->R())*sqr(Tc_)/Pc_;
    	scalar b = 0.07780*this->R()*Tc_/Pc_;
    	scalar Tr = T/Tc_;

    	scalar alpha =
        	sqr
        	(
          	1.0
          	+ (0.37464 + 1.54226*omega_- 0.26992*sqr(omega_))
          	* (1.0 - sqrt(Tr))
        	);

    	scalar B = b*p/(this->R()*T);
    	scalar A = a*alpha*p/sqr(this->R()*T);

    	scalar a2 = B - 1.0;
    	scalar a1 = A - 2.0*B - 3.0*sqr(B);
    	scalar a0 = -A*B + sqr(B) + pow3(B);

    	scalar Q = (3.0*a1 - a2*a2)/9.0;
    	scalar Rl = (9.0*a2*a1 - 27.0*a0 - 2.0*a2*a2*a2)/54;

    	scalar Q3 = Q*Q*Q;
    	scalar D = Q3 + Rl*Rl;

    	scalar root = -1.0;
	
    	if (D <= 0.0)
    	{
        	scalar th = ::acos(Rl/sqrt(-Q3));
        	scalar qm = 2.0*sqrt(-Q);
       	 	scalar r1 = qm*cos(th/3.0) - a2/3.0;
        	scalar r2 = qm*cos((th + 2.0*constant::mathematical::pi)/3.0) - a2/3.0;
        	scalar r3 = qm*cos((th + 4.0*constant::mathematical::pi)/3.0) - a2/3.0;

        	root = max(r1, max(r2, r3));
		rhogas = p/(root*this->R()*T);
			
    	}
    	else
   	{
        	//one root is real
        	scalar D05 = sqrt(D);
        	scalar S = pow(Rl + D05, 1.0/3.0);
        	scalar Tl = 0;
        	if (D05 > Rl)
       	 	{
            		Tl = -pow(mag(Rl - D05), 1.0/3.0);
        	}
        	else
        	{
            		Tl = pow(Rl - D05, 1.0/3.0);
        	}
        
		root = S + Tl - a2/3.0;
    		rhogas = p/(root*this->R()*T);
	}
    }*/
	        
    {
    	double Tgas,pgas,e=NULL,h=NULL,s=NULL,Cv=NULL,Cp=NULL,w=NULL,hjt=NULL;
    	double xgas[]={1.0};
    	long kph = -2, kguess = 0;
    	

	// enthalpy iteration by Newton's method may fail
	if (T < 0.0 || T > 350.0) Info << "abnormal behaviour in enthalpy iteration ...................." << nl << endl;
    
	Tgas = T;
	rhogas = this->rho(p,T)/this->W();
 
	//TPRHOdll(Tgas,pgas,xgas,kph,kguess,rhogas,ierr,herr,255);
    	//THERMdll(Tgas,rhogas,xgas,pgas,e,h,s,Cv,Cp,w,hjt);
	ENTHALdll(Tgas,rhogas,xgas,h);

	hs_ = h*1000.0;

	//Info << "density = " << rhogas << nl << endl;
	//Info << "temperature = " << T << nl << endl;
	//Info << "pressure = " << p << nl << endl;
	//Info << "R = " << this->R() << nl << endl;    	
	//Info << "hs = " << (h*1000.0/this->W()) << nl << endl;	
    	//Info << "Cp*T = " << Cp_*T << nl << endl;
    
	// saturated properties only 
	//hs_ = (-0.000508*Foam::pow (Tgas,4)+0.402133*Foam::pow(Tgas,3)-121.645905*Foam::pow(Tgas,2)+17114.230257*Tgas-537903.887493)*this->W();
	
	//Info << "iterating enthalpy .............." << hs_ << nl << endl;
	//Info << "iterating temperature .............." << Tgas << nl << endl;
	//Info << "iterting pressure .............." << p << nl << endl;

    	//return Cp_*T;
    }

    return (hs_);	
}


template<class equationOfState>
inline Foam::scalar Foam::hRefpropGasThermo<equationOfState>::hc() const
{
    return Hf_;
}


template<class equationOfState>
inline Foam::scalar Foam::hRefpropGasThermo<equationOfState>::s
(
    const scalar p, const scalar T
) const
{

    // Load the shared library
    //-------------------------------------------------------------------
    
    //static int flag = 0;
    //long ierr = 0;
    //char herr[255];

    //if (flag == 0)
    //{
    	//std::string err;
    	//bool loaded_REFPROP = load_REFPROP(err);
    	//long nc = 1, n_PREOS = 2; 
    	//char hfld[] = "CO2.FLD", hhmx[] = "HMX.BNC", href[] = "DEF";

    	//SETUPdll(nc,hfld,hhmx,href,ierr,herr,10000,255,3,255);
        
	// use PR-EoS
	//PREOSdll(n_PREOS);
    	//if (ierr != 0)
    	//{
    		//Info<< "error massage:" << herr << nl << endl;  
   	//}
    	
	//flag++;
    //}
    //-------------------------------------------------------------------		
    
    //double Tgas,rhogas,pgas,e=NULL,h=NULL,s=NULL,cv=NULL,cp=NULL,w=NULL,hjt=NULL,dummy=NULL;
    //double xgas[]={1.0};

    //Tgas = T;
    //rhogas = rho/(this->W());
    //THERM0dll(Tgas,rhogas,xgas,pgas,e,h,s,cv,cp,w,hjt,dummy);

    notImplemented
    (
        "scalar hRefpropGasThermo<equationOfState>::"
        "s(const scalar p, const scalar T) const"
    );
    return T;
    //return (s*1000.0); 
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline void Foam::hRefpropGasThermo<equationOfState>::operator+=
(
    const hRefpropGasThermo<equationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    equationOfState::operator+=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    Cp_ = molr1*Cp_ + molr2*ct.Cp_;
    Hf_ = molr1*Hf_ + molr2*ct.Hf_;
}


template<class equationOfState>
inline void Foam::hRefpropGasThermo<equationOfState>::operator-=
(
    const hRefpropGasThermo<equationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    equationOfState::operator-=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    Cp_ = molr1*Cp_ - molr2*ct.Cp_;
    Hf_ = molr1*Hf_ - molr2*ct.Hf_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::hRefpropGasThermo<equationOfState> Foam::operator+
(
    const hRefpropGasThermo<equationOfState>& ct1,
    const hRefpropGasThermo<equationOfState>& ct2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
      + static_cast<const equationOfState&>(ct2)
    );

    return hRefpropGasThermo<equationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.Cp_
      + ct2.nMoles()/eofs.nMoles()*ct2.Cp_,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
      + ct2.nMoles()/eofs.nMoles()*ct2.Hf_
    );
}


template<class equationOfState>
inline Foam::hRefpropGasThermo<equationOfState> Foam::operator-
(
    const hRefpropGasThermo<equationOfState>& ct1,
    const hRefpropGasThermo<equationOfState>& ct2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
      - static_cast<const equationOfState&>(ct2)
    );

    return hRefpropGasThermo<equationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.Cp_
      - ct2.nMoles()/eofs.nMoles()*ct2.Cp_,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
      - ct2.nMoles()/eofs.nMoles()*ct2.Hf_
    );
}


template<class equationOfState>
inline Foam::hRefpropGasThermo<equationOfState> Foam::operator*
(
    const scalar s,
    const hRefpropGasThermo<equationOfState>& ct
)
{
    return hRefpropGasThermo<equationOfState>
    (
        s*static_cast<const equationOfState&>(ct),
        ct.Cp_,
        ct.Hf_
    );
}


template<class equationOfState>
inline Foam::hRefpropGasThermo<equationOfState> Foam::operator==
(
    const hRefpropGasThermo<equationOfState>& ct1,
    const hRefpropGasThermo<equationOfState>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
