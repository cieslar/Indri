#pragma once

#include<Definitions.hpp>
#include<cmath>
#include<cstdlib>
#include <limits>
#include<Debug.hpp>

template <typename E> class TECartesianExpression 
{
public:
	double operator[](size_t i) const { return static_cast<E const&>(*this)[i];     }

	// The following overload conversions to E, the template argument type;
	// e.g., for VecExpression<VecSum>, this is a conversion to VecSum.
	operator E&()             { return static_cast<      E&>(*this); }
	operator E const&() const { return static_cast<const E&>(*this); }
};

//static double gfRealTypeNaN = std::numeric_limits<RealType>::quiet_NaN();

class TECartesian : public TECartesianExpression<TECartesian> {
public:
	RealType fX;
	RealType fY;
	RealType fZ;

	TECartesian()
	{
		fX=fY=fZ=0.;
	}

	TECartesian(const RealType X, const RealType Y, const RealType Z)
	{
		fX=X;
		fY=Y;
		fZ=Z;
	}

	double& operator[](size_t i)       
	{ 
		if (i==0) return fX;
		if (i==1) return fY;
		if (i==2) return fZ;
		//ERROR("Requesting index of bounds");
		return fX;
	}

	double  operator[](size_t i) const 
	{
		if (i==0) return fX;
		if (i==1) return fY;
		if (i==2) return fZ;
		//ERROR("Requesting index of bounds");
		return fX;
	}


	template <typename E> TECartesian(TECartesianExpression<E> const& Cart)
	{
		for(size_t i(3);i--;) (*this)[i] = Cart[i];
	}

	template <typename E> TECartesian operator=(TECartesianExpression<E> const& Cart)
	{
		for(size_t i(3);i--;) (*this)[i] = Cart[i];
		return *this;
        }
};

template <typename E1, typename E2> class TECartesianSum : public TECartesianExpression<TECartesianSum<E1, E2> > 
{
protected:
	E1 const& _u;
	E2 const& _v;

public:
	TECartesianSum(TECartesianExpression<E1> const& u, TECartesianExpression<E2> const& v) : _u(u), _v(v) {};
	double operator[](size_t i) const { return _u[i] + _v[i]; }
};
 

template <typename E1, typename E2> TECartesianSum<E1,E2> const operator+(TECartesianExpression<E1> const& u, TECartesianExpression<E2> const& v) 
{
   return TECartesianSum<E1, E2>(u, v);
}


template <typename E1, typename E2> class TECartesianSub : public TECartesianExpression<TECartesianSub<E1, E2> > 
{
protected:
	E1 const& _u;
	E2 const& _v;

public:
	TECartesianSub(TECartesianExpression<E1> const& u, TECartesianExpression<E2> const& v) : _u(u), _v(v) {};
	double operator[](size_t i) const { return _u[i] - _v[i]; }
};
 

template <typename E1, typename E2> TECartesianSub<E1,E2> const operator-(TECartesianExpression<E1> const& u, TECartesianExpression<E2> const& v) 
{
   return TECartesianSub<E1, E2>(u, v);
}



template <typename E1, typename E2> class TECartesianMulLeft : public TECartesianExpression<TECartesianMulLeft<E1, E2> > 
{
protected:
	E1 const& _mul;///<multiplication element
	E2 const& _cart;///<cartesian element

public:
	TECartesianMulLeft(E1 const& u, TECartesianExpression<E2> const& v) : _mul(u), _cart(v) {};
	double operator[](size_t i) const { return _mul*_cart[i]; }
};

template <typename E1, typename E2> TECartesianMulLeft<E1,E2> const operator*(E1 const& u, TECartesianExpression<E2> const& v) 
{
   return TECartesianMulLeft<E1, E2>(u, v);
}


template <typename E1, typename E2> class TECartesianMulRight : public TECartesianExpression<TECartesianMulRight<E1, E2> > 
{
protected:
	E1 const& _cart;///<cartesian element
	E2 const& _mul;///<multiplication element

public:
	TECartesianMulRight(TECartesianExpression<E1> const& v, E2 const& u) : _mul(u), _cart(v) {};
	double operator[](size_t i) const { return _mul*_cart[i]; }
};

template <typename E1, typename E2> TECartesianMulRight<E1,E2> const operator*(TECartesianExpression<E1> const& v, E2 const& u) 
{
   return TECartesianMulRight<E1, E2>(v, u);
}

template <typename E> class TECartesianInv : public TECartesianExpression<TECartesianInv<E> > 
{
protected:
	E const& _u;

public:
	TECartesianInv(TECartesianExpression<E> const& u) : _u(u) {};
	double operator[](size_t i) const { return -_u[i]; }
};
 

template <typename E> TECartesianInv<E> const operator-(TECartesianExpression<E> const& u)
{
   return TECartesianInv<E>(u);
}





struct SCartesian 
{
	RealType fX;
	RealType fY;
	RealType fZ;

	SCartesian(){ fX=fY=fZ=0.; };
	SCartesian(const RealType X, const RealType Y, const RealType Z)
	{
		fX=X;
		fY=Y;
		fZ=Z;
	}

};

SCartesian operator+(const SCartesian &c1, const SCartesian &c2);
SCartesian operator-(const SCartesian &c1, const SCartesian &c2);
SCartesian operator-(const SCartesian &c1);
SCartesian operator*(const RealType &c1, const SCartesian &c2);
SCartesian operator*(const SCartesian &c2, const RealType &c1);



//typedef SCartesian Cartesian;
typedef TECartesian Cartesian;




struct Cylindrical
{
	RealType fRho;
	RealType fTheta;
	RealType fZ;
	Cylindrical() { fRho=fTheta=fZ=0.; }
};

struct Galactic
{
	RealType fL;
	RealType fB;
	RealType fDist;
	Galactic() { fL=fB=fDist=0.; }
};

struct Equatorial
{
	RealType fRA;
	RealType fDC;
	RealType fDist;
	Equatorial() { fRA=fDC=fDist=0.; }
};

RealType Arg2PI(const RealType fX, const RealType fY);
Equatorial Gal2Equ(const Galactic Gal);

class Coordinates
{
public:
	Cartesian*   Cart;
	Cylindrical* Cyli;

	Coordinates();
	~Coordinates();

	void Zero();
//	void Cyli2Cart();
	void RecomputeAll();
};


Cartesian Cyli2Cart(const Cylindrical Cyli);
Cartesian RotateVectorXYPlane(const Cartesian Vec, const RealType fAngle);

Cartesian NormVector(const Cartesian Vec);

RealType Length(const Cartesian Vec);

RealType DotProduct(const Cartesian A, const Cartesian B);

Cartesian CrossProduct(const Cartesian A, const Cartesian B);

//Do sprawdzania czy kat l jest w prawo (+1.) czy w lewo (-1.)
RealType GalCheckSide(const Cartesian Vec);

///fSunAngle - angle between ey(0,1,0) and the SunPos vector. [rad]
Galactic Cart2Gal(const Cartesian Pos, const Cartesian SunPos, const RealType fSunAngle);
Galactic Cart2Gal(const Cartesian Pos, const Cartesian SunPos);
