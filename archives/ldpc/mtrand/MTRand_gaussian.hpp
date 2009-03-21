#ifndef INCLUDED_MTRAND_GAUSSIAN_HPP
#define INCLUDED_MTRAND_GAUSSIAN_HPP

#include "mtrand.hpp"

// Gaussian random number generator based on Mersenne twister
class MTRand_gaussian : public MTRand_closed
{
public:
	// Default constructor
	MTRand_gaussian()
	{
		set = false;
	}
	// Constructor with single seed
	MTRand_gaussian(unsigned long seed) : MTRand_closed(seed)
	{
		set = false;
	}
	// Constructor with seed array
	MTRand_gaussian(const unsigned long* seed, int size) : MTRand_closed(seed, size)
	{
		set = false;
	}
	// Generate a random number
	double operator()();

private:
	MTRand_gaussian(const MTRand_gaussian&); // copy constructor not defined
	void operator=(const MTRand_gaussian&); // assignment operator not defined

	bool set;
	double gset;
};

#endif
