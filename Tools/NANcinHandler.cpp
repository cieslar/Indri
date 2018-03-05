#include<NANcinHandler.hpp>

#include<sstream> 
#include<cmath>
#include<limits>
#include<algorithm>


RealType Check4Nan(const string Input)
{
	if(Input.find("nan") != string::npos)
	{
		return numeric_limits<double>::quiet_NaN();
	}

	stringstream Converter;
	Converter<<Input;

	RealType fResult;
	Converter>>fResult;

	return fResult;
}

RealType EmitNan()
{
	return numeric_limits<double>::quiet_NaN();


}

bool IsPresent(const string Input)
{
	if(Input.find("nan") != string::npos)
	{
		return false;
	}

	return true;
}

