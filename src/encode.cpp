/*
* $URL$
* $Date$
* $Rev$
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "ldpcstate.hpp"

using namespace std;

// The SNRs to try.
const double snrs[] =
{
	SNRMIN+0.0,
	SNRMIN+0.1,
	SNRMIN+0.2,
	SNRMIN+0.3,
	SNRMIN+0.4,
	SNRMIN+0.5,
	SNRMIN+0.6,
	SNRMIN+0.7,
	SNRMIN+0.8,
	SNRMIN+0.9
};

///////////////////////////////////////////////////////////////////////////////
// Functors ///////////////////////////////////////////////////////////////////

// Checks the product of the Preaching matrix and the generated parity bits to
// verify them.
struct functor_multhpp {

	bool (&msprod)[M*Z];
	bool success;	// Whether or not the encoding succeeded.

	functor_multhpp(bool (&msprodInit)[M*Z]) : msprod(msprodInit) {}
	
	inline void callbackProduct(int y, bool p) {
		success &= msprod[y] == p;	// Checks if message sum product = parity SP
	}
};

// Functor to find the sum over elements in ms for a row, iterating in x.
struct functor_summsy {
	bool (&ms)[K*Z];
	bool sum;			// The sum for this row.

	functor_summsy(bool (&msInit)[K*Z]) : ms(msInit) {}

	inline void callbackY(int y, int x) {
		sum ^= ms[x];		// Binary addition is equivalent to XOR.
	}
};

///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

LDPCstate::LDPCstate() :
//	grand((unsigned long)time(0)),
//	irand((unsigned long)time(0)),
	mr(H),
	mq(H),
	ms((bool(&)[K*Z])*mx),
	mp((bool(&)[M*Z])mx[K*Z])
{
}

void LDPCstate::init(DecodeMethod::Enum methodInit,
	int snrindexInit)
{
	method = methodInit;
	snrindex = snrindexInit;
	snrdb = snrs[snrindex];
	snr = pow((double)10.0, (double)snrdb/10);
	sigma = pow((double)2.0*RATE*snr, (double)-0.5);
}

void LDPCstate::calculateRho()
{
	// Calculate rho factors for Preaching matrix and submatrices
	int hx = 0, hy = 0,
		hsx = 0, hsy = 0,
		hpx = 0, hpy = 0;
	for (int y = 0; y < M; y++) {
		int hsyi = 0;
		for (int x = 0; x < K; x++)
			if (Hs.H[y][x] != -1) hsyi++;
		if (hsy < hsyi) hsy = hsyi;
	}
	for (int x = 0; x < K; x++) {
		int hsxi = 0;
		for (int y = 0; y < M; y++)
			if (Hs.H[y][x] != -1) hsxi++;
		if (hsx < hsxi) hsx = hsxi;
	}
	for (int a = 0; a < M; a++) {
		int hpxi = 0, hpyi = 0;
		for (int b = 0; b < M; b++) {
			if (Hp.H[a][b] != -1) hpyi++;
			if (Hp.H[b][a] != -1) hpxi++;
		}
		if (hpx < hpxi) hpx = hpxi;
		if (hpy < hpyi) hpy = hpyi;
	}
	for (int y = 0; y < M; y++) {
		int hyi = 0;
		for (int x = 0; x < N; x++)
			if (H.H[y][x] != -1) hyi++;
		if (hy < hyi) hy = hyi;
	}
	for (int x = 0; x < N; x++) {
		int hxi = 0;
		for (int y = 0; y < M; y++)
			if (H.H[y][x] != -1) hxi++;
		if (hx < hxi) hx = hxi;
	}
	
	cout << "M="  << M   << " N="   << N   << " Variant=" << VARIANT << '\n'
		<< "hx="  << hx  << " hy="  << hy  << '\n'
		<< "hsx=" << hsx << " hsy=" << hsy << '\n'
		<< "hpx=" << hpx << " hpy=" << hpy << '\n';
}

void LDPCstate::encode()
{
	// Generate the random bits of the message
	for (int m = 0; m < K*Z; m++)
		ms[m] = irand() & 1;

	// Get the parity bits
	setParity();

#ifdef _DEBUG
	// Double-check that the encoding succeeded
	Hs.multCol(ms, msprod);

	functor_multhpp func(msprod);
	func.success = true;
	Hp.multColCallback(func, mp);

	if (!func.success)
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
void LDPCstate::setParity()
{
	// Magic number for left side of v(0) determination
	// With our matrix, this element is ZERO therefore there is NO SHIFT NEEDED
	// const int xshift = 0;

	functor_summsy funcsum(ms);

	// Determine v(0)
	for (int mi = 0; mi < Z; mi++) // Iterate over the index of v0
	{
		funcsum.sum = 0;
		for (int m = mi; m < Z*M; m += Z)	// Iterate over m for whole H matrix
			Hs.iterY(funcsum, m);			// Effectively iterate over n
		mp[mi] = funcsum.sum;
	}

	// Determine v(1)
	for (int mi = 0; mi < Z; mi++)
	{
		funcsum.sum = Hp.pshift(0,0,mp,mi);	// P(i,k)v(0)
		Hs.iterY(funcsum, mi);				// sigma P(i,j)u(j)
		mp[Z+mi] = funcsum.sum;				// p(1)
	}

	// Determine v(i)
	bool *pmp = mp+Z;	// p(i) starting at i=1
	for (int i = 1; i <= M-2; i++)
	{
		for (int mi = 0; mi < Z; mi++, pmp++)
		{
			funcsum.sum = *pmp ^ Hp.pshift(i,0,mp,mi);	// v(i) + P(i,k)v(0)
			Hs.iterY(funcsum, Z*i+mi);					// sigma P(i,j)u(j)
			pmp[Z] = funcsum.sum;						// p(i+1)
		}
	}
}
