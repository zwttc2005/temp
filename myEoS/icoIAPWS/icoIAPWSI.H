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

#include "icoIAPWS.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::icoIAPWS<Specie, PolySize>::icoIAPWS
(
    const Specie& sp,
    const Polynomial<PolySize>& rhoCoeffs
)
:
    Specie(sp),
    rhoCoeffs_(rhoCoeffs)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::icoIAPWS<Specie, PolySize>::icoIAPWS
(
    const icoIAPWS<Specie, PolySize>& ip
)
:
    Specie(ip),
    rhoCoeffs_(ip.rhoCoeffs_)
{}


template<class Specie, int PolySize>
inline Foam::icoIAPWS<Specie, PolySize>::icoIAPWS
(
    const word& name,
    const icoIAPWS<Specie, PolySize>& ip
)
:
    Specie(name, ip),
    rhoCoeffs_(ip.rhoCoeffs_)
{}


template<class Specie, int PolySize>
inline Foam::autoPtr<Foam::icoIAPWS<Specie, PolySize> >
Foam::icoIAPWS<Specie, PolySize>::clone() const
{
    return autoPtr<icoIAPWS<Specie, PolySize> >
    (
        new icoIAPWS<Specie, PolySize>(*this)
    );
}


template<class Specie, int PolySize>
inline Foam::autoPtr<Foam::icoIAPWS<Specie, PolySize> >
Foam::icoIAPWS<Specie, PolySize>::New(Istream& is)
{
    return autoPtr<icoIAPWS<Specie, PolySize> >
    (
        new icoIAPWS<Specie, PolySize>(is)
    );
}


template<class Specie, int PolySize>
inline Foam::autoPtr<Foam::icoIAPWS<Specie, PolySize> >
Foam::icoIAPWS<Specie, PolySize>::New(const dictionary& dict)
{
    return autoPtr<icoIAPWS<Specie, PolySize> >
    (
        new icoIAPWS<Specie, PolySize>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::scalar Foam::icoIAPWS<Specie, PolySize>::rho
(
    scalar,
    scalar T
) const
{
    return 0.0; // TODO: needs changing
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoIAPWS<Specie, PolySize>::psi
(
    scalar,
    scalar
) const
{
    return 0.0; // TODO: needs changing
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoIAPWS<Specie, PolySize>::Z
(
    scalar,
    scalar
) const
{
    return 0.0; // TODO: needs changing
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoIAPWS<Specie, PolySize>::cpMcv
(
    scalar p,
    scalar T
) const
{
    return 0.0; // TODO: needs changing
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::icoIAPWS<Specie, PolySize>&
Foam::icoIAPWS<Specie, PolySize>::operator=
(
    const icoIAPWS<Specie, PolySize>& ip
)
{
    Specie::operator=(ip);

    rhoCoeffs_ = ip.rhoCoeffs_;

    return *this;
}


template<class Specie, int PolySize>
inline void Foam::icoIAPWS<Specie, PolySize>::operator+=
(
    const icoIAPWS<Specie, PolySize>& ip
)
{
    scalar molr1 = this->nMoles();

    Specie::operator+=(ip);

    molr1 /= this->nMoles();
    scalar molr2 = ip.nMoles()/this->nMoles();

    rhoCoeffs_ = molr1*rhoCoeffs_ + molr2*ip.rhoCoeffs_;
}


template<class Specie, int PolySize>
inline void Foam::icoIAPWS<Specie, PolySize>::operator-=
(
    const icoIAPWS<Specie, PolySize>& ip
)
{
    scalar molr1 = this->nMoles();

    Specie::operator-=(ip);

    molr1 /= this->nMoles();
    scalar molr2 = ip.nMoles()/this->nMoles();

    rhoCoeffs_ = molr1*rhoCoeffs_ - molr2*ip.rhoCoeffs_;
}


template<class Specie, int PolySize>
inline void Foam::icoIAPWS<Specie, PolySize>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
Foam::icoIAPWS<Specie, PolySize> Foam::operator+
(
    const icoIAPWS<Specie, PolySize>& ip1,
    const icoIAPWS<Specie, PolySize>& ip2
)
{
    scalar nMoles = ip1.nMoles() + ip2.nMoles();
    scalar molr1 = ip1.nMoles()/nMoles;
    scalar molr2 = ip2.nMoles()/nMoles;

    return icoIAPWS<Specie, PolySize>
    (
        static_cast<const Specie&>(ip1)
      + static_cast<const Specie&>(ip2),
        molr1*ip1.rhoCoeffs_ + molr2*ip2.rhoCoeffs_
    );
}


template<class Specie, int PolySize>
Foam::icoIAPWS<Specie, PolySize> Foam::operator-
(
    const icoIAPWS<Specie, PolySize>& ip1,
    const icoIAPWS<Specie, PolySize>& ip2
)
{
    scalar nMoles = ip1.nMoles() + ip2.nMoles();
    scalar molr1 = ip1.nMoles()/nMoles;
    scalar molr2 = ip2.nMoles()/nMoles;

    return icoIAPWS<Specie, PolySize>
    (
        static_cast<const Specie&>(ip1)
      - static_cast<const Specie&>(ip2),
        molr1*ip1.rhoCoeffs_ - molr2*ip2.rhoCoeffs_
    );
}


template<class Specie, int PolySize>
Foam::icoIAPWS<Specie, PolySize> Foam::operator*
(
    const scalar s,
    const icoIAPWS<Specie, PolySize>& ip
)
{
    return icoIAPWS<Specie, PolySize>
    (
        s*static_cast<const Specie&>(ip),
        ip.rhoCoeffs_
    );
}


template<class Specie, int PolySize>
Foam::icoIAPWS<Specie, PolySize> Foam::operator==
(
    const icoIAPWS<Specie, PolySize>& ip1,
    const icoIAPWS<Specie, PolySize>& ip2
)
{
    return ip2 - ip1;
}


// ************************************************************************* //
