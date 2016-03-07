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
inline Foam::hRefpropLiqThermo<equationOfState>::hRefpropLiqThermo
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
inline Foam::hRefpropLiqThermo<equationOfState>::hRefpropLiqThermo
(
    const word& name,
    const hRefpropLiqThermo& ct
)
:
    equationOfState(name, ct),
    Cp_(ct.Cp_),
    Hf_(ct.Hf_)
{}


template<class equationOfState>
inline Foam::autoPtr<Foam::hRefpropLiqThermo<equationOfState> >
Foam::hRefpropLiqThermo<equationOfState>::clone() const
{
    return autoPtr<hRefpropLiqThermo<equationOfState> >
    (
        new hRefpropLiqThermo<equationOfState>(*this)
    );
}


template<class equationOfState>
inline Foam::autoPtr<Foam::hRefpropLiqThermo<equationOfState> >
Foam::hRefpropLiqThermo<equationOfState>::New(Istream& is)
{
    return autoPtr<hRefpropLiqThermo<equationOfState> >
    (
        new hRefpropLiqThermo<equationOfState>(is)
    );
}


template<class equationOfState>
inline Foam::autoPtr<Foam::hRefpropLiqThermo<equationOfState> >
Foam::hRefpropLiqThermo<equationOfState>::New(const dictionary& dict)
{
    return autoPtr<hRefpropLiqThermo<equationOfState> >
    (
        new hRefpropLiqThermo<equationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class EquationOfState>
inline Foam::scalar Foam::hRefpropLiqThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class equationOfState>
inline Foam::scalar Foam::hRefpropLiqThermo<equationOfState>::cp
(
    const scalar p,
    const scalar T

) const
{    
    scalar cp_;	 
    scalar rholiq;

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

        	root = min(r1, min(r2, r3));
		rholiq = p/(root*this->R()*T);
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
		rholiq = p/(root*this->R()*T);
    	}
    }*/
    
    {
    	double Tliq,pliq,e=NULL,h=NULL,s=NULL,Cv=NULL,Cp=NULL,w=NULL,hjt=NULL;
    	double xliq[]={1.0};
    	//long kph = -1, kguess = 0; // integer required for TP-calculation

    	Tliq = T;        //K
	rholiq = this->rho(p,T)/this->W();
    	//pliq = p/1000.0; //kpa

        
	//TPRHOdll(Tliq,pliq,xliq,kph,kguess,rholiq,ierr,herr,255);	
	
    	//THERMdll(Tliq,rholiq,xliq,pliq,e,h,s,Cv,Cp,w,hjt);
	
	CVCPdll(Tliq,rholiq,xliq,Cv,Cp);

    	//if (p < 2500000)
    	//{
      		//Info << "internal density = " << rholiq << nl << endl;
    		//Info << "internal temperature = " << Tliq << nl << endl;
      		//Info << "hs = " << (h*1000.0/44.01) << nl << endl;  
    		//Info << "Cp*T = " << Cp_*T << nl << endl;
    	//}
	
	cp_ = Cp*1000.0;
    }
    
    return (cp_);
}

template<class equationOfState>
inline Foam::scalar Foam::hRefpropLiqThermo<equationOfState>::ha
(
    const scalar p, const scalar T
) const
{
    return Cp_*T + Hf_;
}


template<class equationOfState>
inline Foam::scalar Foam::hRefpropLiqThermo<equationOfState>::hs
(
    const scalar p, const scalar T
) const
{
    scalar hs_;
    scalar rholiq;	

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

        	root = min(r1, min(r2, r3));
		rholiq = p/(root*this->R()*T);
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
		rholiq = p/(root*this->R()*T);
    	}
    }*/    
    {
    	double Tliq,pliq,e=NULL,h=NULL,s=NULL,Cv=NULL,Cp=NULL,w=NULL,hjt=NULL;
    	double xliq[]={1.0};
    	//long kph = -1, kguess = 0;
    
    	Tliq = T;
	rholiq = this->rho(p,T)/this->W();
    	
	//TPRHOdll(Tliq,pliq,xliq,kph,kguess,rholiq,ierr,herr,255);
    	//THERMdll(Tliq,rholiq,xliq,pliq,e,h,s,Cv,Cp,w,hjt);	
	ENTHALdll(Tliq,rholiq,xliq,h);

	hs_ = h*1000.0;
    	
	//Info << "hs = " << (h*1000.0/this->W()) << nl << endl;	
	//Info << "liquid density = " << rholiq*this->W() << nl << endl;	
	//Info << "liquid T = " << Tliq << nl << endl;    	
	//Info << "Cp*T = " << Cp_*T << nl << endl;
    
	// saturated properties only 
	//hs_ = (0.000491*Foam::pow(Tliq,4)-0.397733*Foam::pow(Tliq,3)+120.103770*Foam::pow(Tliq,2)-14045.518453*Tliq+449986.580984)*this->W();


    	//return Cp_*T;
    }

    return (hs_);	
}


template<class equationOfState>
inline Foam::scalar Foam::hRefpropLiqThermo<equationOfState>::hc() const
{
    return Hf_;
}


template<class equationOfState>
inline Foam::scalar Foam::hRefpropLiqThermo<equationOfState>::s
(
    const scalar rho, const scalar T
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
    
    //double Tliq,rholiq,p = NULL,e=NULL,h=NULL,s=NULL,cv=NULL,cp=NULL,w=NULL,hjt=NULL,dummy=NULL;
    //double xliq[]={1.0};

    //Tliq = T;
    //rholiq = rho/(this->W());
    //THERM0dll(Tliq,rholiq,xliq,p,e,h,s,cv,cp,w,hjt,dummy);

    notImplemented
    (
        "scalar hRefpropLiqThermo<equationOfState>::"
        "s(const scalar p, const scalar T) const"
    );
    return T;
    //return (s*1000.0); 
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline void Foam::hRefpropLiqThermo<equationOfState>::operator+=
(
    const hRefpropLiqThermo<equationOfState>& ct
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
inline void Foam::hRefpropLiqThermo<equationOfState>::operator-=
(
    const hRefpropLiqThermo<equationOfState>& ct
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
inline Foam::hRefpropLiqThermo<equationOfState> Foam::operator+
(
    const hRefpropLiqThermo<equationOfState>& ct1,
    const hRefpropLiqThermo<equationOfState>& ct2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
      + static_cast<const equationOfState&>(ct2)
    );

    return hRefpropLiqThermo<equationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.Cp_
      + ct2.nMoles()/eofs.nMoles()*ct2.Cp_,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
      + ct2.nMoles()/eofs.nMoles()*ct2.Hf_
    );
}


template<class equationOfState>
inline Foam::hRefpropLiqThermo<equationOfState> Foam::operator-
(
    const hRefpropLiqThermo<equationOfState>& ct1,
    const hRefpropLiqThermo<equationOfState>& ct2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
      - static_cast<const equationOfState&>(ct2)
    );

    return hRefpropLiqThermo<equationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.Cp_
      - ct2.nMoles()/eofs.nMoles()*ct2.Cp_,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
      - ct2.nMoles()/eofs.nMoles()*ct2.Hf_
    );
}


template<class equationOfState>
inline Foam::hRefpropLiqThermo<equationOfState> Foam::operator*
(
    const scalar s,
    const hRefpropLiqThermo<equationOfState>& ct
)
{
    return hRefpropLiqThermo<equationOfState>
    (
        s*static_cast<const equationOfState&>(ct),
        ct.Cp_,
        ct.Hf_
    );
}


template<class equationOfState>
inline Foam::hRefpropLiqThermo<equationOfState> Foam::operator==
(
    const hRefpropLiqThermo<equationOfState>& ct1,
    const hRefpropLiqThermo<equationOfState>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
