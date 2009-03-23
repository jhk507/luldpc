/*
* $URL$
* $Date$
* $Rev$

Warning: This code is VERY NON-REENTRANT. This is deliberate to facilitate
access to variables from functors and to greatly increase performance by
cutting down on frame pointer generation. In short, no multithreading allowed.
*/

#include <ctime>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>

#define OUTPUT_DEBUGFILE 0	// Enable to output data to a debug file

#include "mtrand/MTRand_gaussian.hpp"
#include "ldpc.hpp"

using namespace std;

namespace LDPC
{

// Data for unexpanded half-rate Preaching matrix H
const int Ha[M][N] =
{
	{-1,94,73,-1,-1,-1,-1,-1,55,83,-1,-1, 7, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{-1,27,-1,-1,-1,22,79, 9,-1,-1,-1,12,-1, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{-1,-1,-1,24,22,81,-1,33,-1,-1,-1, 0,-1,-1, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1},
	{61,-1,47,-1,-1,-1,-1,-1,65,25,-1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,-1,-1},
	{-1,-1,39,-1,-1,-1,84,-1,-1,41,72,-1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,-1},
	{-1,-1,-1,-1,46,40,-1,82,-1,-1,-1,79, 0,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1},
	{-1,-1,95,53,-1,-1,-1,-1,-1,14,18,-1,-1,-1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1},
	{-1,11,73,-1,-1,-1, 2,-1,-1,47,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,-1,-1,-1},
	{12,-1,-1,-1,83,24,-1,43,-1,-1,-1,51,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,-1,-1},
	{-1,-1,-1,-1,-1,94,-1,59,-1,-1,70,72,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,-1},
	{-1,-1, 7,65,-1,-1,-1,-1,39,49,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0},
	{43,-1,-1,-1,-1,66,-1,41,-1,-1,-1,26, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0}
};

// Matrix sparsity parameters
#define RHO_H_Y  7
#define RHO_H_X  6
#define RHO_HS_Y 5
#define RHO_HS_X 6
#define RHO_HP_Y 3
#define RHO_HP_X 3

const Preaching<M,N,RHO_H_Y, RHO_H_X>	H(Ha, 0);	// Unexpanded half-rate Preaching matrix H
const Preaching<M,K,RHO_HS_Y,RHO_HS_X>	Hs(Ha, 0);	// Unexpanded half-rate Preaching matrix H (first half)
const Preaching<M,M,RHO_HP_Y,RHO_HP_X>	Hp(Ha, K);	// Unexpanded half-rate Preaching matrix H (second half, for parity)

Automatrix1<bool, N*Z> mx;			// (col) Combination of ms and mp
Automatrix1<long double, N*Z> my;	// (col) Encoder output after AWGN

// Set the aliases into mx
bool (&ms)[K*Z] = (bool(&)[K*Z])*mx.getData(0);		// (col) Message
bool (&mp)[M*Z] = (bool(&)[M*Z])*mx.getData(K*Z);	// (col) Generated parity

Automatrix1<bool, M*Z> msprod;		// Encoding verification column

PreachingBased<long double, M,N,RHO_H_Y,RHO_H_X> mr(H);		// R matrix
PreachingBased<long double, M,N,RHO_H_Y,RHO_H_X> mq(H);		// Q matrix
PreachingBased<long double, M,N,RHO_H_Y,RHO_H_X> mq0(H);	// Q matrix (iteration 0)
Automatrix1<long double, N*Z> ml;	// L column
Automatrix1<long double, N*Z> ml0;	// L column (iteration 0)
Automatrix1<bool, N*Z> mxhat;		// xhat column

// The Gaussian distribution random number generator
MTRand_gaussian grand((unsigned long)time(0));
// Discrete value random number generator
MTRand_int32 irand((unsigned long)~time(0));

// Constants for AWGN calculation
const double R = 0.5;	// Rate
double snr;				// Signal-to-noise ratio
double snrdb;			// Signal-to-noise ratio (decibels)
double sigma;

// Beta (for minsum decoding)
const double beta = 0.15;

#if OUTPUT_DEBUGFILE
ofstream debugfile("debugfile.tsv");
#endif

int imax;

const enum DecodeMethod
{
	bp,
	offms
} method = offms;

void init()
{
/*	// Calculate rho factors for Preaching matrix and submatrices
	int hx = 0, hy = 0,
		hsx = 0, hsy = 0,
		hpx = 0, hpy = 0;
	for (int a = 0; a < M; a++) {
		int hsxi = 0, hsyi = 0,
			hpxi = 0, hpyi = 0;
		for (int b = 0; b < M; b++) {
			if (Hs.H[a][b] != -1) hsyi++;
			if (Hs.H[b][a] != -1) hsxi++;
			if (Hp.H[a][b] != -1) hpyi++;
			if (Hp.H[b][a] != -1) hpxi++;
		}
		if (hsx < hsxi) hsx = hsxi;
		if (hsy < hsyi) hsy = hsyi;
		if (hpx < hpxi) hpx = hpxi;
		if (hpy < hpyi) hpy = hpyi;
	}
	for (int a = 0; a < M; a++) {
		int hyi = 0;
		for (int b = 0; b < N; b++)
			if (H.H[a][b] != -1) hyi++;
		if (hy < hyi) hy = hyi;
	}
	for (int a = 0; a < N; a++) {
		int hxi = 0;
		for (int b = 0; b < M; b++)
			if (H.H[b][a] != -1) hxi++;
		if (hx < hxi) hx = hxi;
	}
*/

	cout << "Enter signal to noise ratio (dB): ";
//	cin >> snrdb;
	snrdb = 1.5;
	snr = pow(10.0, snrdb/10);
	sigma = pow(2*R*snr, -0.5);

	cout << "Enter maximum decoding iteration count: ";
//	cin >> imax;
	imax = 50;
}

void execute()
{
	// Initialize the simulation
	init();

	int nerrs = 0;	// The number of block errors

#if OUTPUT_DEBUGFILE
	debugfile << "Unexpanded half-rate Preaching matrix" << endl;
	for (int m = 0; m < M; m++)
	{
		for (int n = 0; n < N; n++)
			debugfile << Ha[m][n] << '\t';
		debugfile << endl;
	}
	debugfile << endl;

	debugfile << "Expanded half-rate Preaching matrix" << endl;
	H.output(debugfile);
	
	debugfile.flush();
#endif

	// The main block loop
	for (int b = 1; nerrs < 50; b++)
	{
		// Encode
		encode();

#if OUTPUT_DEBUGFILE
		debugfile << "Message:" << endl;
		outputLargeContiguous<bool,K,Z>(ms, debugfile);

		debugfile << "Encoded parity bits:" << endl;
		outputLargeContiguous<bool,M,Z>(mp, debugfile);
#endif

		// Decode
		if (!decode())
			nerrs++;

		cout << "Block " << b << ": " << nerrs << " errors, BLER=" << 100.0*nerrs/b << '%' << endl;
#if OUTPUT_DEBUGFILE
		debugfile << "Block " << b << ": " << nerrs << " errors, BLER=" << 100.0*nerrs/b << '%' << endl;
		break;
#endif
	}
}


void encode()
{
	// Generate the random bits of the message
	for (int m = 0; m < K*Z; m++)
		ms[m] = irand() & 1;

	// Encode
#ifdef _DEBUG
	cout << "Encoding..." << endl;
#if OUTPUT_DEBUGFILE
	debugfile << "Encoding..." << endl;
#endif
#endif

	// Get the parity bits
	setParity();

#ifdef _DEBUG
	// Double-check that the encoding succeeded
	Hs.multCol(ms, msprod);

	static bool success;
	success = true;

	struct functor_multhpp {
		static inline void callback(int y, bool p) {
			success &= msprod[y] == p;	// Checks if message sum product = parity SP
		}
	};
	Hp.multCol<functor_multhpp>(mp);
	
	cout << "Encoder check " << (success ? "passed." : "failed.") << endl;
#if OUTPUT_DEBUGFILE
	debugfile << "Encoder check " << (success ? "passed." : "failed.") << endl;
#endif
#endif

	for (int n = 0; n < N*Z; n++)
		// Perform BPSK and AWGN addition
		my[n] = (mx[n] ? -1 : 1) + grand()*sigma;	// 1->-1 and 0->1
}

// Set the parity matrix based on the message
void setParity()
{
	static bool sum;

	// Functor to find the sum over elements in ms for a row,
	// iterating in x
	struct functor_summsy {
		static inline void callbackY(int y, int x, int xi) {
			sum ^= ms[x]; //XOR
		}
	};

	// Magic number for left side of v(0) determination
	// With our matrix, this element is ZERO therefore there is NO SHIFT NEEDED
	// const int xshift = 0;

	// Determine v(0)
	for (int mi = 0; mi < Z; mi++) // Iterate over the index of v0
	{
		sum = 0;
		for (int m = mi; m < Z*M; m += Z)	// Iterate over m for whole H matrix
			Hs.iterY<functor_summsy>(m);	// Effectively iterate over n
		mp[mi] = sum;
	}

	// Determine v(1)
	for (int mi = 0; mi < Z; mi++)
	{
		sum = Hp.pshift(0,0,mp,mi);		// P(i,k)v(0)
		Hs.iterY<functor_summsy>(mi);	// sigma P(i,j)u(j)
		mp[Z+mi] = sum;					// p(1)
	}

	// Determine v(i)
	bool *pmp = mp+Z;	// p(i) starting at i=1
	for (int i = 1; i <= M-2; i++)
	{
		for (int mi = 0; mi < Z; mi++, pmp++)
		{
			sum = *pmp ^ Hp.pshift(i,0,mp,mi);	// v(i) + P(i,k)v(0)
			Hs.iterY<functor_summsy>(Z*i+mi);	// sigma P(i,j)u(j)
			pmp[Z] = sum;						// p(i+1)
		}
	}
}


bool decode()
{
	// Set initial state
	for (int n = 0; n < N*Z; n++)
	{
		static long double l;

		if (method == bp)
			l = 2.0/sigma/sigma*my[n]; // Required for BP algorithm
		else
			l = my[n];

		ml0[n] = l;
		ml[n] = l;

		struct functor_setq {
			static inline void callbackX(int y, int x, int yi) {
				mq0.Vxc[x][yi] = l;
				mq.Vxc[x][yi] = l;
			}
		};
		H.iterX<functor_setq>(n);
	}

	// Iterative decoding
	for (int i = 0; ; )
	{
#if OUTPUT_DEBUGFILE
		debugfile << "Before iteration " << i << ":" << endl;

		debugfile << "l:" << endl;
		ml.outputLarge<N,Z>(debugfile);

		debugfile << "q:" << endl;
		mq.output(debugfile);

		debugfile.flush();

		if (i >= 1)
			return true;
#endif

		// Update mr
		switch (method)
		{
		case bp:
			rupdate_bp();
			break;
		case offms:
			rupdate_offms();
			break;
		}

#if OUTPUT_DEBUGFILE
		debugfile << "After iteration " << i << ":" << endl;

		debugfile << "r:" << endl;
		mr.output(debugfile);

		debugfile.flush();
#endif

		// Update mq and ml
		for (int n = 0; n < Z*N; n++)
		{
			static long double sigma;
			sigma = 0;

			struct functor_sigmar {
				static inline void callbackX(int y, int x, int yi) {
					sigma += mr.Vxc[x][yi];
				}
			};
			H.iterX<functor_sigmar>(n);

			struct functor_updateq
			{
				static inline void callbackX(int y, int x, int yi) {
					// Performs exclusion
					mq.Vxc[x][yi] = mq0.Vxc[x][yi] + sigma - mr.Vxc[x][yi];
				}
			};

			ml[n] = ml0[n] + sigma;

			mxhat[n] = ml[n] < 0; // Hard decision
		}

		// Check that the decoding succeeded
		static int nerrs;	// The number of H*xhat errors
		nerrs = 0;

		struct functor_multhxhat {
			static inline void callback(int y, bool p) {
				nerrs += p;
			}
		};
		H.multCol<functor_multhxhat>(mxhat);

		int diff = 0;	 // The number of x==xhat errors
		for (int j = 0; j < Z*K; j++)
			diff += mxhat[j] != ms[j];

#ifdef _DEBUG
		cout << "Iteration " << setw(3) << i << ": "
			<< setw(4) << nerrs << " H*xhat errors, "
			<< setw(4) << diff << " x==xhat errors."
			<< endl;
#if OUTPUT_DEBUGFILE
		debugfile << "Iteration " << setw(3) << i << ": "
			<< setw(4) << nerrs << " H*xhat errors, "
			<< setw(4) << diff << " x==xhat errors."
			<< endl;
#endif
#endif

		if (!nerrs)
			return true;

		if (++i > imax)
		{
			cout << "Maximum iteration reached; error." << endl;
#if OUTPUT_DEBUGFILE
			debugfile << "Maximum iteration reached; error." << endl;
#endif
			return false;
		}
	}
}

void rupdate_bp()
{
	// Update mr
	for (int m = 0; m < Z*M; m++)
	{
		// Do the graph iteration to calculate the pi term without exclusion
		static long double pi;
		pi = 1;

		struct functor_r_bp_pi {
			static inline void callbackY(int y, int x, int xi) {
				pi *= tanh(mq.Vyc[y][xi]/2);
			}
		};
		H.iterY<functor_r_bp_pi>(m);

		struct functor_r_bp_update {
			static inline void callbackY(int y, int x, int xi) {
				long double pir = pi;
				const long double tanhr = tanh(mq.Vyc[y][xi]/2);
				if (tanhr)
					pir /= tanhr;
				else
					exit(-1);	// Divide by 0
				if (pir == 1)
				{
//					mr[m][n] = numeric_limits<long double>::max();
					exit(-1);	// Divide by 0
				}
				else
				{
					const long double lnarg = (1+pir)/(1-pir);
					if (lnarg <= 0)
						exit(-1);	// Negative log
					else
						mr.Vyc[y][xi] = log(lnarg);
				}
			}
		};
		H.iterY<functor_r_bp_update>(m);
	}
}

void rupdate_offms()
{
	// Update mr
	// OFF-MS method
	for (int m = 0; m < Z*M; m++)
	{
		static int pi;
		pi = 1;	// Multiplicative identity

		// Min function identities
		// Assume a very large value so that it may be overwritten on the first
		// iteration.
		static long double min0, min1;
		min0 = numeric_limits<long double>::max();
		min1 = min0;

		// Do the graph iteration to calculate the pi and min
		// terms without exclusion
		struct functor_r_offms_pi {
			static inline void callbackY(int y, int x, int xi) {
				long double qv = mq.Vyc[y][xi];
				if (qv < 0)
					pi = -pi;
				qv = fabs(qv);
				if (min0 >= qv)
				{
					min1 = min0;
					min0 = qv;
				}
				else if (min1 > qv)
					min1 = qv;
			}
		};
		H.iterY<functor_r_offms_pi>(m);

		struct functor_r_offms_update {
			static inline void callbackY(int y, int x, int xi) {
				int pir = pi;

				const long double qv = mq.Vyc[y][xi];
				// Perform exclusion on the pi term
				if (qv < 0)
					pir = -pir;
				// Perform exclusion on the min term
				const long double qvmin = (fabs(qv) == min0) ? min1 : min0;

				// Offset min sum calculation for r^(i)_(m,n)
				mr.Vyc[y][xi] = pir * max((long double)0.0, qvmin - beta);
			}
		};
		H.iterY<functor_r_offms_update>(m);
	}
}

}
