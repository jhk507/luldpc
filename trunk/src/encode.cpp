/*
* $URL$
* $Date$
* $Rev$
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "mtrand/MTRand_gaussian.hpp"

#include "encode.hpp"

using namespace std;

namespace LDPC
{

///////////////////////////////////////////////////////////////////////////////
// Globals ////////////////////////////////////////////////////////////////////

// The Gaussian distribution random number generator
MTRand_gaussian grand(0);	//((unsigned long)time(0));
// Discrete value random number generator
MTRand_int32 irand(1);		//((unsigned long)~time(0));

bool mx[N*Z];	// (col) Combination of ms and mp
double my[N*Z];	// (col) Encoder output after AWGN

// Set the aliases into mx
bool (&ms)[K*Z] = (bool(&)[K*Z])mx;			// (col) Message
bool (&mp)[M*Z] = (bool(&)[M*Z])mx[K*Z];	// (col) Generated parity

bool msprod[M*Z];	// Encoding verification column

// The SNRs to try.
const double snrs[] =
{
	1.0,
	1.1,
	1.2,
	1.3,
	1.4,
	1.5,
	1.6,
	1.7,
	1.8,
	1.9
};

int snrindex;								// The current SNR index

// Constants for AWGN calculation
double snr;	// Signal-to-noise ratio
double sigma;

///////////////////////////////////////////////////////////////////////////////
// Functors ///////////////////////////////////////////////////////////////////

// Checks the product of the Preaching matrix and the generated parity bits to
// verify them.
struct functor_multhpp {
	static bool success;	// Whether or not the encoding succeeded.
	static inline void callbackProduct(int y, bool p) {
		success &= msprod[y] == p;	// Checks if message sum product = parity SP
	}
};
bool functor_multhpp::success;

// Functor to find the sum over elements in ms for a row, iterating in x.
struct functor_summsy {
	static bool sum;		// The sum for this row.
	static inline void callbackY(int y, int x) {
		sum ^= ms[x];		// Binary addition is equivalent to XOR.
	}
};
bool functor_summsy::sum;

///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

void setSnrDB(double snrdb)
{
	cout << "Setting the SNR to " << snrdb << " dB...\n";
	// Calculate the linear SNR and sigma from an SNR in dB
	snr = pow((double)10.0, (double)snrdb/10);
	sigma = pow((double)2.0*RATE*snr, (double)-0.5);
}

void encode()
{
	// Generate the random bits of the message
	for (int m = 0; m < K*Z; m++)
		ms[m] = irand() & 1;

	// Get the parity bits
	setParity();

#ifdef _DEBUG
	// Double-check that the encoding succeeded
	Hs.multCol(ms, msprod);

	functor_multhpp::success = true;
	Hp.multCol<functor_multhpp>(mp);

	if (!functor_multhpp::success)
	{
		cerr << "Encoding check failed!\n";
		exit(-1);
	}
#endif

	for (int n = 0; n < N*Z; n++)
		// Perform BPSK and AWGN addition
		my[n] = (mx[n] ? -1 : 1) + grand()*sigma;	// 1->-1 and 0->1
}

// Set the parity matrix based on the message
void setParity()
{
	// Magic number for left side of v(0) determination
	// With our matrix, this element is ZERO therefore there is NO SHIFT NEEDED
	// const int xshift = 0;

	// Determine v(0)
	for (int mi = 0; mi < Z; mi++) // Iterate over the index of v0
	{
		functor_summsy::sum = 0;
		for (int m = mi; m < Z*M; m += Z)	// Iterate over m for whole H matrix
			Hs.iterY<functor_summsy>(m);	// Effectively iterate over n
		mp[mi] = functor_summsy::sum;
	}

	// Determine v(1)
	for (int mi = 0; mi < Z; mi++)
	{
		functor_summsy::sum = Hp.pshift(0,0,mp,mi);	// P(i,k)v(0)
		Hs.iterY<functor_summsy>(mi);				// sigma P(i,j)u(j)
		mp[Z+mi] = functor_summsy::sum;				// p(1)
	}

	// Determine v(i)
	bool *pmp = mp+Z;	// p(i) starting at i=1
	for (int i = 1; i <= M-2; i++)
	{
		for (int mi = 0; mi < Z; mi++, pmp++)
		{
			functor_summsy::sum = *pmp ^ Hp.pshift(i,0,mp,mi);	// v(i) + P(i,k)v(0)
			Hs.iterY<functor_summsy>(Z*i+mi);					// sigma P(i,j)u(j)
			pmp[Z] = functor_summsy::sum;						// p(i+1)
		}
	}
}

}
