/*
* $URL: https://luldpc.googlecode.com/svn/trunk/src/preachingbased.hpp $
* $Date: 2009-03-25 21:01:53 -0400 (Wed, 25 Mar 2009) $
* $Rev: 30 $
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
