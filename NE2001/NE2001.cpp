#include<NE2001.hpp>


float fort_dmdsm(float l, float b, int ndir, float dmpsr, float dist)
{
	float sm,smtau,smtheta,smiso;
	char limit;
	int limit_size=1;
	FORTRAN_NAME(dmdsm)(&l, &b, &ndir, &dmpsr, &dist, &limit, &limit_size, &sm, &smtau, &smtheta, &smiso);

	if(ndir >= 0) 	return dist;
	else		return dmpsr;
}

