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
#include <cstdlib>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>

#define OUTPUT_DEBUGFILE 0	// Enable to output data to a debug file

#include "mtrand/MTRand_gaussian.hpp"
#include "histogram.hpp"
#include "preachingbased.hpp"
#include "ldpc.hpp"

using namespace std;

namespace LDPC
{

///////////////////////////////////////////////////////////////////////////////
// Globals ////////////////////////////////////////////////////////////////////

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

bool mx[N*Z];			// (col) Combination of ms and mp
long double my[N*Z];	// (col) Encoder output after AWGN

// Set the aliases into mx
bool (&ms)[K*Z] = (bool(&)[K*Z])mx;			// (col) Message
bool (&mp)[M*Z] = (bool(&)[M*Z])mx[K*Z];	// (col) Generated parity

bool msprod[M*Z];		// Encoding verification column

// Decoding matrices
PreachingBased<long double, M,N,RHO_H_Y,RHO_H_X> mr(H);		// R matrix
PreachingBased<long double, M,N,RHO_H_Y,RHO_H_X> mq(H);		// Q matrix
long double ml[N*Z];	// L column
long double ml0[N*Z];	// L column (iteration 0)
bool mxhat[N*Z];		// xhat column

// The maximum number of decode iterations
const int imax = 50;

// The orthagonality error and output difference error histograms.
// The template parameters are the number of histogram buckets, the full size
// of the data range, and the desired portion of the data range to examine.
Histogram<20, M*Z, (int)(M*Z*0.33)> orthhist[imax];
Histogram<20, N*Z, (int)(N*Z*0.06)> diffhist[imax];

// The Gaussian distribution random number generator
MTRand_gaussian grand(0);	//((unsigned long)time(0));
// Discrete value random number generator
MTRand_int32 irand(1);		//((unsigned long)~time(0));

// Constants for AWGN calculation
long double snr;	// Signal-to-noise ratio
long double sigma;

// Beta (for minsum decoding)
const long double beta = 0.15;

#if OUTPUT_DEBUGFILE
ofstream debugfile("debugfile.tsv");
#endif

// The decode method.
const enum DecodeMethod
{
	bp,		// Belief propagation
	offms	// Offset min sum
} method = offms;

///////////////////////////////////////////////////////////////////////////////
// Functors ///////////////////////////////////////////////////////////////////

struct functor_multhpp {
	static bool success;
	static inline void callbackProduct(int y, bool p) {
		success &= msprod[y] == p;	// Checks if message sum product = parity SP
	}
};
bool functor_multhpp::success;

// Functor to find the sum over elements in ms for a row,
// iterating in x
struct functor_summsy {
	static bool sum;
	static inline void callbackY(int y, int x) {
		sum ^= ms[x]; //XOR
	}
};
bool functor_summsy::sum;

// Check that the decoding succeeded with orthagonality verification

// The matrix multiplication functor
struct functor_multhxhat {
	static int nerrs;
	static inline void callbackProduct(int y, bool p) {
		// There is an error every time there is a 1 in the product.
		nerrs += p;
	}
};
int functor_multhxhat::nerrs;	// The number of H*xhat errors

// The functor to set the initial values for Q and Q0
struct functor_setq {
	static long double l;
	static inline void callback(long double &q) {
		q = l;
	}
};
long double functor_setq::l;


// Do the graph iteration to calculate the pi term without exclusion
struct functor_r_bp_pi {
	static long double pi;					// The pi term (without exclusion)
	static long double tanh_cache[RHO_H_Y];	// The cache of calculated tanh values
	static int xi;							// The current index in the row

	static inline void callback(long double &q) {
		const long double tanhval = tanh(q/2.0);
		pi *= tanhval;
		tanh_cache[xi++] = tanhval;
	}
};
long double functor_r_bp_pi::pi;
long double functor_r_bp_pi::tanh_cache[RHO_H_Y];
int functor_r_bp_pi::xi;

// The functor to update the R matrix
struct functor_r_bp_update {
	static inline void callback(long double &r, long double &q) {
		long double pir = functor_r_bp_pi::pi;
		const long double tanhr = functor_r_bp_pi::tanh_cache[functor_r_bp_pi::xi++];
		if (tanhr)
			pir /= tanhr;
		else
		{
			pir = numeric_limits<long double>::max();
			cerr << "Warning: Divide by 0 in BP for pir!\n";
//					exit(-1);	// Divide by 0
		}
		if (pir == 1)
		{
			r = numeric_limits<long double>::max();
			cerr << "Warning: Divide by 0 in BP for r!\n";
//					exit(-1);	// Divide by 0
		}
		else
		{
			const long double lnarg = (1+pir)/(1-pir);
			if (lnarg <= 0)
			{
				r = -numeric_limits<long double>::max();
				cerr << "Warning: Negative log in BP!";
//						exit(-1);	// Negative log
			}
			else
				r = log(lnarg);
		}
	}
};

// Do the graph iteration to calculate the pi and min
// terms without exclusion
struct functor_r_offms_pi {
	static long double min0, min1;
	static int pi;
	static inline void callback(long double &q) {
		long double qv = q;
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
long double functor_r_offms_pi::min0;
long double functor_r_offms_pi::min1;
int functor_r_offms_pi::pi;

struct functor_r_offms_update {
	static inline void callback(long double &r, long double &q) {
		int pir = functor_r_offms_pi::pi;

		const long double qv = q;
		// Perform exclusion on the pi term
		if (qv < 0)
			pir = -pir;
		// Perform exclusion on the min term
		const long double qvmin = (fabs(qv) == functor_r_offms_pi::min0) ?
				functor_r_offms_pi::min1 : functor_r_offms_pi::min0;

		// Offset min sum calculation for r^(i)_(m,n)
		r = pir * max((long double)0.0, qvmin - beta);
	}
};

struct functor_sigmar {
	static long double rsigma;
	static inline void callback(long double &r) {
		// Calculates the sigma term without exclusion
		rsigma += r;
	}
};
long double functor_sigmar::rsigma;

struct functor_updateq {
	static long double q0;
	static inline void callback(long double &q, long double &r) {
		// Performs exclusion and sets Q
		q = q0 + functor_sigmar::rsigma - r;
	}
};
long double functor_updateq::q0;


///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

void calculateRho()
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
}

void setSnrDB(long double snrdb)
{
	// Calculate the linear SNR and sigma from an SNR in dB
	snr = pow((long double)10.0, (long double)snrdb/10);
	sigma = pow((long double)2.0*RATE*snr, (long double)-0.5);
}

void execute()
{
	// Initialize the simulation
	setSnrDB(1.5);

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
	for (int b = 1; b <= 500; b++)
	{
		// Encode
		encode();

#if OUTPUT_DEBUGFILE
		debugfile << "Message:" << endl;
		outputLargeContiguous<K,Z>(ms, debugfile);

		debugfile << "Encoded parity bits:" << endl;
		outputLargeContiguous<M,Z>(mp, debugfile);
#endif

		// Decode
		if (!decode())
			nerrs++;

		if (!(b%10))
		{
			cout << "Block errors: " << nerrs << " / " << b << "\tBLER=" << 100.0*nerrs/b << '%' << endl;
#if OUTPUT_DEBUGFILE
			debugfile << "Block errors: " << nerrs << " / " << b << "\tBLER=" << 100.0*nerrs/b << '%' << endl;
			break;
#endif
		}
	}

	// Output the error histograms.
	ofstream hist("histogram.tsv");

	// i is on the vertical axis, buckets are on the horizontal axis.
	hist << '\t';
	orthhist->outputHeader(hist);
	for (int i = 0; i < imax; i++)
	{
		hist << i << '\t';
		orthhist[i].output(hist);
	}

	hist << "\n\n\t";
	diffhist->outputHeader(hist);
	for (int i = 0; i < imax; i++)
	{
		hist << i << '\t';
		diffhist[i].output(hist);
	}
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
		functor_summsy::sum = Hp.pshift(0,0,mp,mi);		// P(i,k)v(0)
		Hs.iterY<functor_summsy>(mi);	// sigma P(i,j)u(j)
		mp[Z+mi] = functor_summsy::sum;					// p(1)
	}

	// Determine v(i)
	bool *pmp = mp+Z;	// p(i) starting at i=1
	for (int i = 1; i <= M-2; i++)
	{
		for (int mi = 0; mi < Z; mi++, pmp++)
		{
			functor_summsy::sum = *pmp ^ Hp.pshift(i,0,mp,mi);	// v(i) + P(i,k)v(0)
			Hs.iterY<functor_summsy>(Z*i+mi);	// sigma P(i,j)u(j)
			pmp[Z] = functor_summsy::sum;						// p(i+1)
		}
	}
}


bool decode()
{
	// Set the initial state of the decoder
	decode_initial();

	// Iterative decoding
	for (int i = 0; ; )
	{
#if OUTPUT_DEBUGFILE
		debugfile << "Before iteration " << i << ":" << endl;

		debugfile << "l:" << endl;
		outputLarge<N,Z>(ml, debugfile);

		debugfile << "q:" << endl;
		mq.output(debugfile);

		debugfile.flush();

		if (i >= 1)
			return true;
#endif

		// Update the R matrix
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

		// Update the Q and L matrices
		qlupdate();

		functor_multhxhat::nerrs = 0;
		H.multCol<functor_multhxhat>(mxhat);
		orthhist[i].report(functor_multhxhat::nerrs);

		int diff = 0;	 // The number of x==xhat errors
		for (int j = 0; j < Z*K; j++)
			diff += mxhat[j] != ms[j];
		diffhist[i].report(diff);

		if (!functor_multhxhat::nerrs)
		{
			if (diff)
			{
				cerr << "Warning: False positive; " << diff << " errors.\n";
				return false;
			}
			for (++i; i < imax; i++)
			{
				orthhist[i].report(0);
				diffhist[i].report(0);
			}
			return true;
		}

		if (++i > imax)
			return false;
	}
}

void decode_initial()
{
	// Set initial state
	for (int n = 0; n < N*Z; n++)
	{

		if (method == bp)
			functor_setq::l = 2.0/sigma/sigma*my[n]; // Required for BP algorithm
		else
			functor_setq::l = my[n];

		ml0[n] = functor_setq::l;
		ml[n] = functor_setq::l;

		mq.iterX<functor_setq>(n);
	}
}

void rupdate_bp()
{
	// Update mr
	for (int m = 0; m < Z*M; m++)
	{
		functor_r_bp_pi::pi = 1;
		functor_r_bp_pi::xi = 0;
		mq.iterY<functor_r_bp_pi>(m);

		// Reset the row index
		functor_r_bp_pi::xi = 0;
		mr.iterY2<functor_r_bp_update>(m,mq);
	}
}

void rupdate_offms()
{
	// Update mr
	// OFF-MS method
	for (int m = 0; m < Z*M; m++)
	{
		functor_r_offms_pi::pi = 1;	// Multiplicative identity

		// Min function identities
		// Assume a very large value so that it may be overwritten on the first
		// iteration.
		functor_r_offms_pi::min0 = numeric_limits<long double>::max();
		functor_r_offms_pi::min1 = numeric_limits<long double>::max();

		mq.iterY<functor_r_offms_pi>(m);

		mr.iterY2<functor_r_offms_update>(m,mq);
	}
}

void qlupdate()
{
	// Update mq and ml
	for (int n = 0; n < Z*N; n++)
	{
		functor_sigmar::rsigma = 0;

		mr.iterX<functor_sigmar>(n);

		functor_updateq::q0 = ml0[n];

		mq.iterX2<functor_updateq>(n,mr);

		ml[n] = ml0[n] + functor_sigmar::rsigma;

		mxhat[n] = ml[n] < 0; // Hard decision
	}
}

}
