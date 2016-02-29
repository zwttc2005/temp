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
    const scalar rho,
    const scalar T

) const
{
    
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
    
    double Tgas,rhogas,p = NULL,e=NULL,h=NULL,s=NULL,cv=NULL,cp=NULL,w=NULL,hjt=NULL,dummy=NULL;
    double xgas[]={1.0};

    Tgas = T;
    rhogas = rho/(this->W());
    THERM0dll(Tgas,rhogas,xgas,p,e,h,s,cv,cp,w,hjt,dummy);

    
    //Info << "hs = " << (h*1000.0) << nl << endl;  
    //Info << "Cp*T = " << Cp_*T << nl << endl;
    return (cp*1000.0);
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
    const scalar rho, const scalar T
) const
{
    
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
    
    double Tgas,rhogas,p = NULL,e=NULL,h=NULL,s=NULL,cv=NULL,cp=NULL,w=NULL,hjt=NULL,dummy=NULL;
    double xgas[]={1.0};

    Tgas = T;
    rhogas = rho/(this->W());
    THERM0dll(Tgas,rhogas,xgas,p,e,h,s,cv,cp,w,hjt,dummy);

    
    //Info << "hs = " << (h*1000.0) << nl << endl;  
    //Info << "Cp*T = " << Cp_*T << nl << endl;
    
    //return Cp_*T;
    return (h*1000.0);	
}


template<class equationOfState>
inline Foam::scalar Foam::hRefpropGasThermo<equationOfState>::hc() const
{
    return Hf_;
}


template<class equationOfState>
inline Foam::scalar Foam::hRefpropGasThermo<equationOfState>::s
(
    const scalar rho, const scalar T
) const
{

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
    
    double Tgas,rhogas,p = NULL,e=NULL,h=NULL,s=NULL,cv=NULL,cp=NULL,w=NULL,hjt=NULL,dummy=NULL;
    double xgas[]={1.0};

    Tgas = T;
    rhogas = rho/(this->W());
    THERM0dll(Tgas,rhogas,xgas,p,e,h,s,cv,cp,w,hjt,dummy);

    //notImplemented
    //(
        //"scalar hRefpropGasThermo<equationOfState>::"
        //"s(const scalar p, const scalar T) const"
    //);
    //return T;
    return (s*1000.0); 
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
