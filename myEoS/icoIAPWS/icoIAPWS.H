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
    Foam::icoIAPWS

Description
    Incompressible, polynomial form of equation of state, using a polynomial
    function for density based on IAPWS-IF97.

    WARNING: This is just a dummy class, the functionallity still needs to be
    implemented!

SourceFiles
    icoIAPWSI.H
    icoIAPWS.C

\*---------------------------------------------------------------------------*/

#ifndef icoIAPWS_H
#define icoIAPWS_H

#include "autoPtr.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie, int PolySize>
class icoIAPWS;

template<class Specie, int PolySize>
icoIAPWS<Specie, PolySize> operator+
(
    const icoIAPWS<Specie, PolySize>&,
    const icoIAPWS<Specie, PolySize>&
);

template<class Specie, int PolySize>
icoIAPWS<Specie, PolySize> operator-
(
    const icoIAPWS<Specie, PolySize>&,
    const icoIAPWS<Specie, PolySize>&
);

template<class Specie, int PolySize>
icoIAPWS<Specie, PolySize> operator*
(
    const scalar,
    const icoIAPWS<Specie, PolySize>&
);

template<class Specie, int PolySize>
icoIAPWS<Specie, PolySize> operator==
(
    const icoIAPWS<Specie, PolySize>&,
    const icoIAPWS<Specie, PolySize>&
);

template<class Specie, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const icoIAPWS<Specie, PolySize>&
);


/*---------------------------------------------------------------------------*\
                        Class icoIAPWS Declaration
\*---------------------------------------------------------------------------*/

template<class Specie, int PolySize=8>
class icoIAPWS
:
    public Specie
{
    // Private data

        //- Density polynomial coefficients
        //  Note: input in [kg/m3], but internally uses [kg/m3/kmol]
        Polynomial<PolySize> rhoCoeffs_;


public:

    // Constructors

        //- Construct from components
        inline icoIAPWS
        (
            const Specie& sp,
            const Polynomial<PolySize>& rhoPoly
        );

        //- Construct from Istream
        icoIAPWS(Istream&);

        //- Construct from dictionary
        icoIAPWS(const dictionary& dict);

        //- Construct as copy
        inline icoIAPWS(const icoIAPWS&);

        //- Construct as named copy
        inline icoIAPWS(const word& name, const icoIAPWS&);

        //- Construct and return a clone
        inline autoPtr<icoIAPWS> clone() const;

        // Selector from Istream
        inline static autoPtr<icoIAPWS> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<icoIAPWS> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "icoIAPWS<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kmol K]
            inline scalar cpMcv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline icoIAPWS& operator=(const icoIAPWS&);
        inline void operator+=(const icoIAPWS&);
        inline void operator-=(const icoIAPWS&);

        inline void operator*=(const scalar);


    // Friend operators

        friend icoIAPWS operator+ <Specie, PolySize>
        (
            const icoIAPWS&,
            const icoIAPWS&
        );

        friend icoIAPWS operator- <Specie, PolySize>
        (
            const icoIAPWS&,
            const icoIAPWS&
        );

        friend icoIAPWS operator* <Specie, PolySize>
        (
            const scalar s,
            const icoIAPWS&
        );

        friend icoIAPWS operator== <Specie, PolySize>
        (
            const icoIAPWS&,
            const icoIAPWS&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie, PolySize>
        (
            Ostream&,
            const icoIAPWS&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeIcoPolynomial(PolySize)                                          \
                                                                             \
defineTemplateTypeNameAndDebugWithName                                       \
(                                                                            \
    icoIAPWS<Specie, PolySize>,                                         \
    "icoIAPWS<"#PolySize">",                                            \
    0                                                                        \
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "icoIAPWSI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "icoIAPWS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
