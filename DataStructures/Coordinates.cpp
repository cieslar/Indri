#include<Coordinates.hpp>
#include<cmath>
#include<Debug.hpp>
#include<iostream>
using namespace std;

SCartesian operator+(const SCartesian &c1, const SCartesian &c2)
{
	SCartesian Result;


	Result.fZ=c1.fZ+c2.fZ;
	Result.fX=c1.fX+c2.fX;
	Result.fY=c1.fY+c2.fY;

    	return Result;
}

SCartesian operator-(const SCartesian &c1, const SCartesian &c2)
{
	SCartesian Result;
	Result.fX=c1.fX-c2.fX;
	Result.fY=c1.fY-c2.fY;
	Result.fZ=c1.fZ-c2.fZ;
    	return Result;
}
SCartesian operator-(const SCartesian &c1)
{
	SCartesian Result;
	Result.fX=-c1.fX;
	Result.fY=-c1.fY;
	Result.fZ=-c1.fZ;
    	return Result;
}
SCartesian operator*(const RealType &c1, const SCartesian &c2)
{
	SCartesian Result;
	Result.fX=c1*c2.fX;
	Result.fY=c1*c2.fY;
	Result.fZ=c1*c2.fZ;
    	return Result;
}
SCartesian operator*(const SCartesian &c2, const RealType &c1)
{
	SCartesian Result;
	Result.fX=c1*c2.fX;
	Result.fY=c1*c2.fY;
	Result.fZ=c1*c2.fZ;
    	return Result;
}




Coordinates::~Coordinates()
{
	delete Cart;
	
	
	delete Cyli;
}

Coordinates::Coordinates()
{
	Cart=new Cartesian;
	Cyli=new Cylindrical;
	Zero();
}

void Coordinates::Zero()
{
	Cart->fX=Cart->fY=Cart->fZ=0.;
	Cyli->fRho=Cyli->fTheta=Cyli->fZ=0.;
}

Cartesian Cyli2Cart(const Cylindrical Cyli)
{
        Cartesian Cart;
        Cart.fX=Cyli.fRho*cos(Cyli.fTheta);
        Cart.fY=Cyli.fRho*sin(Cyli.fTheta);
        Cart.fZ=Cyli.fZ;
        return Cart;
}



Cartesian RotateVectorXYPlane(const Cartesian Vec, const RealType fAngle)
{
	Cartesian Res;
	Res.fX=Vec.fX*cos(fAngle) - Vec.fY*sin(fAngle); 
	Res.fY=Vec.fX*sin(fAngle) + Vec.fY*cos(fAngle);
	Res.fZ=Vec.fZ;
	return Res; 
}

Cartesian NormVector(const Cartesian Vec)
{
	Cartesian Res;
	RealType fLength=sqrt(Vec.fX*Vec.fX + Vec.fY*Vec.fY + Vec.fZ*Vec.fZ);
	Res.fX=Vec.fX/fLength;
	Res.fY=Vec.fY/fLength;
	Res.fZ=Vec.fZ/fLength;
	return Res;

}

RealType Length(const Cartesian Vec)
{
	return sqrt(Vec.fX*Vec.fX + Vec.fY*Vec.fY + Vec.fZ*Vec.fZ);
}

RealType DotProduct(const Cartesian A, const Cartesian B)
{
	return A.fX*B.fX + A.fY*B.fY + A.fZ*B.fZ;
}

Cartesian CrossProduct(const Cartesian A, const Cartesian B)
{
	Cartesian Result;

	Result.fX = A.fY*B.fZ - A.fZ*B.fY;
	Result.fY = A.fZ*B.fX - A.fX*B.fZ;
	Result.fZ = A.fX*B.fY - A.fY*B.fX;

	return Result;
}

//Do sprawdzania czy kat l jest w prawo (+1.) czy w lewo (-1.)
RealType GalCheckSide(const Cartesian Vec)
{
	//zetowa skladowa: Vec x ey(0,1,0)
	if( (Vec.fX*1. - Vec.fY) >= 0. )
	{
		return 1.;
	}
	return -1.;
}
Galactic Cart2Gal(const Cartesian Pos, const Cartesian SunPos)
{
	Galactic Gal;
	Cartesian LocalPos=RotateVectorXYPlane(Pos-SunPos,M_PI/2.);

	RealType fDistPlane = sqrt(LocalPos.fX*LocalPos.fX + LocalPos.fY*LocalPos.fY);

	Gal.fDist=Length(LocalPos);
	Gal.fB = atan2(LocalPos.fZ, fDistPlane);

	if(Gal.fB > M_PI/2. || Gal.fB < -M_PI/2.)
	{
		ERROR("Coordinates' transformation is wrong.");
	}

	Gal.fL = atan2(LocalPos.fY,LocalPos.fX);
	if(Gal.fL<0.) Gal.fL+=2.*M_PI;

	return Gal;
}

#if 0
///fSunAngle - angle between ey(0,1,0) and the SunPos vector. [rad]
Galactic Cart2Gal(const Cartesian Pos, const Cartesian SunPos, const RealType fSunAngle)
{
	Cartesian LocalPos, ey;
	Galactic Gal;	

	ey.fX=ey.fZ=0.;
	ey.fY=1.;

	LocalPos=Pos-SunPos;

	Gal.fDist=Length(LocalPos);



	Gal.fL=acos(DotProduct(LocalPos,ey)/Gal.fDist);
	RealType fZ=LocalPos.fZ;
	LocalPos.fZ=0.;
	Gal.fL*=GalCheckSide(LocalPos);

	Gal.fB=atan(fZ/Length(LocalPos));

	return Gal;
}

Galactic Cart2Gal(const Cartesian Pos, const Cartesian SunPos)
{

	Cartesian LocalPos;
	Galactic Gal;	

	LocalPos = RotateVectorXYPlane(Pos-SunPos, M_PI/2.);


	Gal.fDist=Length(LocalPos);

	RealType fDXY = sqrt( LocalPos.fX*LocalPos.fX + LocalPos.fY*LocalPos.fY);

   	Gal.fL = atan2( (-1.0)*LocalPos.fY,  LocalPos.fX );
	Gal.fB = atan2( LocalPos.fZ, fDXY );


	Gal.fL = Gal.fL - static_cast<int>(Gal.fL/(2.*M_PI)) * 2.*M_PI;


	return Gal;
}
#endif

/* gal2eq.c

   Converts galactic coordinates to equatorial coordinates.
   2000.0

   Formulas and data from the article
   The 1997 reference of diffuse night sky brightness
   by Leinert et al
   A&A Supplement series, Vol. 127, January I 1998, 1-99
   http://aanda.u-strasbg.fr:2002/articles/aas/full/1998/01/ds1449/ds1449.html
   See also: Caroll: 937

   input:
      l = galactic longitude, from ... to, radians
      b 
   output:
      alpha = right ascension angle 0,2pi, radians
      delta - declination angle -pi/2,pi/2, radians

   TO DO
use a define???
    B1950
   l_zero=33*M_PI/180;
   alpha_zero=282.25*M_PI/180;
   delta_NGP=(27.4)*M_PI/180;


*/


Equatorial Gal2Equ(const Galactic Gal)
{

//void gal2eq(const RealType fL, const RealType fB, RealType *alpha, RealType *delta)

   /* conversion constants */
   RealType l_zero;
   RealType alpha_zero;
   RealType delta_NGP;
   /* intermediate values */
   RealType sin_delta, caa, saa, aa;

	RealType alpha, delta;

   /* 2000.0 values */
   l_zero=32.93*M_PI/180;
   alpha_zero=282.86*M_PI/180;
   delta_NGP=(27+7.8/60)*M_PI/180;


   sin_delta = sin(Gal.fB)*sin(delta_NGP) + cos(Gal.fB)*cos(delta_NGP)*sin(Gal.fL-l_zero);
   delta = asin(sin_delta);
   caa = cos(Gal.fL-l_zero)*cos(Gal.fB)/cos(delta);
   saa = ( -sin(Gal.fB)*cos(delta_NGP) + cos(Gal.fB)*sin(delta_NGP)*sin(Gal.fL-l_zero) )/cos(delta);
   aa=Arg2PI(caa, saa);
   alpha=aa+alpha_zero;
   /*
   Because of the previous addition,
   this needs to be reduced again to [0,2pi)
   */
   alpha = alpha - 2*M_PI * floor(alpha/(2*M_PI));

	Equatorial Equ;
	Equ.fRA=alpha;
	Equ.fDC=delta;
	Equ.fDist=Gal.fDist;
	return Equ;
}
int sgn(double x )
{
   return (x>0)-(x<0);
}

RealType Arg2PI(const RealType fX, const RealType fY)
{

   RealType phi;

   /*
   This is here to avoid a useless warning, because gcc tries to be too
   clever (basically trying to run the code at compile time)
   and ends up being dumb.
   */
   phi=0;


   if ( fX == 0 && fY == 0 )
      phi=0;
   if ( fX == 0 && fY != 0 )
      phi=(2- sgn(fY)   )*M_PI/2;
   if ( fX > 0 && fY >= 0 )
      phi=atan(fY/fX);
   if ( fX > 0 && fY < 0 )
      phi=atan(fY/fX)+2*M_PI;
   if ( fX < 0 )
      phi=atan(fY/fX)+M_PI;

   return phi;

}

