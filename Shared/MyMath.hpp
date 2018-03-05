#pragma once

#include<Definitions.hpp>

#define MINIMUM(A,B) (A<B?A:B)
#define MAXIMUM(A,B) (A>B?A:B)


int Factorial(const int n);
RealType LnFactorial(const int n);

RealType LnPoissonProbability(const int nLambda, const int nK);

RealType LnGaussProbability(const RealType fX, const RealType fMu, const RealType fSigma=1.);
