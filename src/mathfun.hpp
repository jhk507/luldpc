/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include <cmath>

namespace LDPC
{

#define LUTSIZE 1024
#define TANHMAX 19.1

extern const double tanhlut[LUTSIZE];

void initmath();

void makeluts();

inline double tanhapp(double x)
{
	const double xabs = fabs(x);
	const double pval = (xabs >= TANHMAX) ? 1 :
		tanhlut[(int)(xabs*(LUTSIZE/TANHMAX))];
	return (x >= 0) ? pval : -pval;
}

};
