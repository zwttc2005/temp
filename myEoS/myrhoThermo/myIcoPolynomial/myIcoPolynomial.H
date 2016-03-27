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
    Foam::myIcoPolynomial

Description
    Incompressible, polynomial form of equation of state, using a polynomial
    function for density.

SourceFiles
    myIcoPolynomialI.H
    myIcoPolynomial.C

\*---------------------------------------------------------------------------*/

#ifndef myIcoPolynomial_H
#define myIcoPolynomial_H

#include "autoPtr.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie, int PolySize>
class myIcoPolynomial;

template<class Specie, int PolySize>
myIcoPolynomial<Specie, PolySize> operator+
(
    const myIcoPolynomial<Specie, PolySize>&,
    const myIcoPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
myIcoPolynomial<Specie, PolySize> operator-
(
    const myIcoPolynomial<Specie, PolySize>&,
    const myIcoPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
myIcoPolynomial<Specie, PolySize> operator*
(
    const scalar,
    const myIcoPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
myIcoPolynomial<Specie, PolySize> operator==
(
    const myIcoPolynomial<Specie, PolySize>&,
    const myIcoPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const myIcoPolynomial<Specie, PolySize>&
);


/*---------------------------------------------------------------------------*\
                        Class myIcoPolynomial Declaration
\*---------------------------------------------------------------------------*/

template<class Specie, int PolySize=8>
class myIcoPolynomial
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
        inline myIcoPolynomial
        (
            const Specie& sp,
            const Polynomial<PolySize>& rhoPoly
        );

        //- Construct from Istream
        myIcoPolynomial(Istream&);

        //- Construct from dictionary
        myIcoPolynomial(const dictionary& dict);

        //- Construct as copy
        inline myIcoPolynomial(const myIcoPolynomial&);

        //- Construct as named copy
        inline myIcoPolynomial(const word& name, const myIcoPolynomial&);

        //- Construct and return a clone
        inline autoPtr<myIcoPolynomial> clone() const;

        // Selector from Istream
        inline static autoPtr<myIcoPolynomial> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<myIcoPolynomial> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "myIcoPolynomial<" + word(Specie::typeName_()) + '>';
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

        inline myIcoPolynomial& operator=(const myIcoPolynomial&);
        inline void operator+=(const myIcoPolynomial&);
        inline void operator-=(const myIcoPolynomial&);

        inline void operator*=(const scalar);


    // Friend operators

        friend myIcoPolynomial operator+ <Specie, PolySize>
        (
            const myIcoPolynomial&,
            const myIcoPolynomial&
        );

        friend myIcoPolynomial operator- <Specie, PolySize>
        (
            const myIcoPolynomial&,
            const myIcoPolynomial&
        );

        friend myIcoPolynomial operator* <Specie, PolySize>
        (
            const scalar s,
            const myIcoPolynomial&
        );

        friend myIcoPolynomial operator== <Specie, PolySize>
        (
            const myIcoPolynomial&,
            const myIcoPolynomial&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie, PolySize>
        (
            Ostream&,
            const myIcoPolynomial&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makemyIcoPolynomial(PolySize)                                          \
                                                                             \
defineTemplateTypeNameAndDebugWithName                                       \
(                                                                            \
    myIcoPolynomial<Specie, PolySize>,                                         \
    "myIcoPolynomial<"#PolySize">",                                            \
    0                                                                        \
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "myIcoPolynomialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "myIcoPolynomial.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
