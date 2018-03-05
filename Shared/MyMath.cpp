#include<MyMath.hpp>
#include<Debug.hpp>
#include<cmath>


int Factorial(const int n)
{
	int nResult(1);
	if(n<0) ERROR("Negative argument for the fatorial.");
	if(n>20) ERROR("Too large argument for the fatorial (over 20).");
	for(int i(2);i<=n;i++) nResult*=i;
	return nResult;
}


RealType LnFactorial(const int n)
{
	if(n<0) ERROR("factorial of negative number!");
	if(n==0) return log(1.);

	RealType fResult(0.);
	for(int i(1);i<=n;i++)
	{
		fResult += log(static_cast<RealType>(i));
	}

	return fResult;
}

#include<iostream>
using namespace std;
RealType LnPoissonProbability(const int nLambda, const int nK)
{
	//cout<<nLambda<<" "<<nK<<endl;
	if(nLambda==0) return -1e10;
	return nK*log(static_cast<double>(nLambda)) - nLambda - LnFactorial(nK);
}

#include<fstream>

RealType LnGaussProbability(const RealType fX, const RealType fMu, const RealType fSigma)
{
	return -log(sqrt(2.*M_PI)*fSigma)-(fX-fMu)*(fX-fMu)/(2.*fSigma*fSigma);
}
