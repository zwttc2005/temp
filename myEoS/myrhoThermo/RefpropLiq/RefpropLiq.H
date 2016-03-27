/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014  OpenFOAM Foundation
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
    Foam::RefpropLiq

Description
    RefpropLiq gas equation of state.

SourceFiles
    RefpropLiqI.H
    RefpropLiq.C

\*---------------------------------------------------------------------------*/

#ifndef RefpropLiq_H
#define RefpropLiq_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class RefpropLiq;

template<class Specie>
inline RefpropLiq<Specie> operator+
(
    const RefpropLiq<Specie>&,
    const RefpropLiq<Specie>&
);

template<class Specie>
inline RefpropLiq<Specie> operator-
(
    const RefpropLiq<Specie>&,
    const RefpropLiq<Specie>&
);

template<class Specie>
inline RefpropLiq<Specie> operator*
(
    const scalar,
    const RefpropLiq<Specie>&
);

template<class Specie>
inline RefpropLiq<Specie> operator==
(
    const RefpropLiq<Specie>&,
    const RefpropLiq<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const RefpropLiq<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class RefpropLiq Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class RefpropLiq
:
    public Specie
{
    // Private data

        //- Critical Temperature [K]
        scalar Tc_;

        //- Critical volume [m^3/kmol]
        scalar Vc_;

        //- Critical compression factor [-]
        scalar Zc_;

        //- Critical Pressure [Pa]
        scalar Pc_;

        //- Acentric factor [-]
        scalar omega_;


public:

    // Constructors

        //- Construct from components
        inline RefpropLiq
        (
            const Specie& sp,
            const scalar& Tc,
            const scalar& Vc,
            const scalar& Zc,
            const scalar& Pc,
            const scalar& omega
        );

        //- Construct from Istream
        RefpropLiq(Istream&);

        //- Construct from dictionary
        RefpropLiq(const dictionary& dict);

        //- Construct as named copy
        inline RefpropLiq(const word& name, const RefpropLiq&);

        //- Construct and return a clone
        inline autoPtr<RefpropLiq> clone() const;

        // Selector from Istream
        inline static autoPtr<RefpropLiq> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<RefpropLiq> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "RefpropLiq<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties


            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor [-]
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kmol K]
            inline scalar cpMcv(scalar p, scalar T) const;


    // IO

            //- Write to Ostream
            void write(Ostream& os) const;

    // Member operators

        inline void operator+=(const RefpropLiq&);
        inline void operator-=(const RefpropLiq&);

        inline void operator*=(const scalar);


    // Friend operators

        friend RefpropLiq operator+ <Specie>
        (
            const RefpropLiq&,
            const RefpropLiq&
        );

        friend RefpropLiq operator- <Specie>
        (
            const RefpropLiq&,
            const RefpropLiq&
        );

        friend RefpropLiq operator* <Specie>
        (
            const scalar s,
            const RefpropLiq&
        );

        friend RefpropLiq operator== <Specie>
        (
            const RefpropLiq&,
            const RefpropLiq&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const RefpropLiq&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "RefpropLiqI.H"

#ifdef NoRepository
#   include "RefpropLiq.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
