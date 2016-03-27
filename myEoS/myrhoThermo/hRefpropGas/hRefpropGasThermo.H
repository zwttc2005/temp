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

Class
    Foam::hRefpropGasThermo

Description
    Constant properties thermodynamics package
    templated into the equationOfState.

SourceFiles
    hRefpropGasI.H
    hRefpropGasThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hRefpropGasThermo_H
#define hRefpropGasThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class hRefpropGasThermo;

template<class EquationOfState>
inline hRefpropGasThermo<EquationOfState> operator+
(
    const hRefpropGasThermo<EquationOfState>&,
    const hRefpropGasThermo<EquationOfState>&
);

template<class EquationOfState>
inline hRefpropGasThermo<EquationOfState> operator-
(
    const hRefpropGasThermo<EquationOfState>&,
    const hRefpropGasThermo<EquationOfState>&
);

template<class EquationOfState>
inline hRefpropGasThermo<EquationOfState> operator*
(
    const scalar,
    const hRefpropGasThermo<EquationOfState>&
);

template<class EquationOfState>
inline hRefpropGasThermo<EquationOfState> operator==
(
    const hRefpropGasThermo<EquationOfState>&,
    const hRefpropGasThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hRefpropGasThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class hRefpropGasThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hRefpropGasThermo
:
    public EquationOfState
{
    // Private data

        scalar Cp_;
        scalar Hf_;


    // Private Member Functions

        //- Construct from components
        inline hRefpropGasThermo
        (
            const EquationOfState& st,
            const scalar cp,
            const scalar hf
        );


public:

    // Constructors

        //- Construct from Istream
        hRefpropGasThermo(Istream&);

        //- Construct from dictionary
        hRefpropGasThermo(const dictionary& dict);

        //- Construct as named copy
        inline hRefpropGasThermo(const word&, const hRefpropGasThermo&);

        //- Construct and return a clone
        inline autoPtr<hRefpropGasThermo> clone() const;

        //- Selector from Istream
        inline static autoPtr<hRefpropGasThermo> New(Istream& is);

        //- Selector from dictionary
        inline static autoPtr<hRefpropGasThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hRefpropGas<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties
            
	    //- Heat capacity at constant pressure [J/(kmol K)]
            inline scalar cp(const scalar p, const scalar T) const;
            //- Absolute Enthalpy [J/kmol]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kmol]
            //inline scalar hs(const scalar p, const scalar T) const;
	    inline scalar hs(const scalar p, const scalar T) const;	

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const hRefpropGasThermo&);
        inline void operator-=(const hRefpropGasThermo&);


    // Friend operators

        friend hRefpropGasThermo operator+ <EquationOfState>
        (
            const hRefpropGasThermo&,
            const hRefpropGasThermo&
        );

        friend hRefpropGasThermo operator- <EquationOfState>
        (
            const hRefpropGasThermo&,
            const hRefpropGasThermo&
        );

        friend hRefpropGasThermo operator* <EquationOfState>
        (
            const scalar,
            const hRefpropGasThermo&
        );

        friend hRefpropGasThermo operator== <EquationOfState>
        (
            const hRefpropGasThermo&,
            const hRefpropGasThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hRefpropGasThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hRefpropGasThermoI.H"

#ifdef NoRepository
#   include "hRefpropGasThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
