#include "MTRand_gaussian.hpp"

#include <cmath>

// Generate a random number
double MTRand_gaussian::operator()() {
	if (set)
	{
		set = false;
		return gset;
	}
	double r, v1, v2;
	do {
		v1 = 2.0 * MTRand_closed::operator()() - 1.0;
		v2 = 2.0 * MTRand_closed::operator()() - 1.0;
		r = v1*v1 + v2*v2;
	} while ((r >= 1.0) || (r == 0.0));
	const double fac = sqrt(-2.0*log(r)/r);
	gset = v1*fac;
	set = true;
	return v2*fac;
}
