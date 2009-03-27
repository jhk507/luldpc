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
#include <iomanip>
#include <fstream>
#include <string>

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

bool mx[N*Z];		// (col) Combination of ms and mp
double my[N*Z];		// (col) Encoder output after AWGN

// Set the aliases into mx
bool (&ms)[K*Z] = (bool(&)[K*Z])mx;			// (col) Message
bool (&mp)[M*Z] = (bool(&)[M*Z])mx[K*Z];	// (col) Generated parity

bool msprod[M*Z];	// Encoding verification column

// Decoding matrices
PreachingBased<double, M,N,RHO_H_Y,RHO_H_X> mr(H);	// R matrix
PreachingBased<double, M,N,RHO_H_Y,RHO_H_X> mq(H);	// Q matrix
double ml[N*Z];		// L column
double ml0[N*Z];	// L column (iteration 0)
bool mxhat[N*Z];	// xhat column

#define IMAX 50		// The maximum number of decode iterations
#define NBUCKETS 25	// The number of histogram buckets
#define NBLOCKS 200	// The number of blocks to run

// The decode method.
const enum DecodeMethod
{
	bp = 0,		// Belief propagation
	offms,		// Offset min sum
	ndecodes	// (The number of decoding algorithms)
} method = offms;

const char *const decodeNames[ndecodes] =
{
	"bp",
	"offms"
};

// The SNRs to try.
const double snrs[] = { 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6 };

#define NSNRS (sizeof(snrs)/sizeof(*snrs))	// The number of SNRs to try
#define DEFAULTSNR 5						// The index of the default SNR.
int snrindex;								// The current SNR index

// The orthagonality error and message error histograms.
// The template parameters are the number of histogram buckets, the full size
// of the data range, and the desired portion of the data range to examine.
Histogram<NBUCKETS, M*Z, (int)(M*Z*0.33), NBLOCKS> orthhist[ndecodes][NSNRS][IMAX];
Histogram<NBUCKETS, N*Z, (int)(N*Z*0.06), NBLOCKS> messhist[ndecodes][NSNRS][IMAX];

// The Gaussian distribution random number generator
MTRand_gaussian grand(0);	//((unsigned long)time(0));
// Discrete value random number generator
MTRand_int32 irand(1);		//((unsigned long)~time(0));

// Constants for AWGN calculation
double snr;	// Signal-to-noise ratio
double sigma;

// Beta (for minsum decoding)
const double beta = 0.15;

#if OUTPUT_DEBUGFILE
ofstream debugfile("debugfile.tsv");
#endif

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

// The matrix multiplication functor for decoding orthagonality verification.
struct functor_multhxhat {
	static int nerrs;		// The number of orthagonality errors.
	static inline void callbackProduct(int y, bool p) {
		// There is an error every time there is a 1 in the product.
		nerrs += p;
	}
};
int functor_multhxhat::nerrs;

// The functor to set the initial values for Q.
struct functor_setq {
	static double l;	// The L-matrix value to take.
	static inline void callback(double &q) {
		q = l;
	}
};
double functor_setq::l;


// Calculate the BP decoding pi term without exclusion. Keeps a cache of tanh
// values as these take a while to get.
struct functor_r_bp_pi {
	static double pi;					// The pi term (without exclusion)
	static double tanh_cache[RHO_H_Y];	// The cache of calculated tanh values
	static int xi;							// The current index in the row

	static inline void callback(double &q) {
		const double tanhval = tanh(q/2.0);	// Calculate the tanh term.
		pi *= tanhval;								// Multiply it into pi.
		tanh_cache[xi++] = tanhval;					// Cache it.
	}
};
double functor_r_bp_pi::pi;
double functor_r_bp_pi::tanh_cache[RHO_H_Y];
int functor_r_bp_pi::xi;

// The functor to update the R matrix, using BP decoding. Performs exclusion.
struct functor_r_bp_update {
	static inline void callback(double &r, double &q) {
		double pir = functor_r_bp_pi::pi;		// The pi term to use.
		// Retrieve the appropriate cached tanh.
		const double tanhr = functor_r_bp_pi::tanh_cache[functor_r_bp_pi::xi++];
		if (tanhr)
			// Perform exclusion.
			pir /= tanhr;
		else
		{
			pir = numeric_limits<double>::max();
#ifdef _DEBUG
			cerr << "Warning: Divide by 0 in BP for pir!\n";
#endif
//			exit(-1);	// Divide by 0
		}
		if (pir == 1)
		{
			r = numeric_limits<double>::max();
#ifdef _DEBUG
			cerr << "Warning: Divide by 0 in BP for r!\n";
#endif
//			exit(-1);	// Divide by 0
		}
		else
		{
			// Calculate the argument for ln.
			const double lnarg = (1+pir)/(1-pir);
			if (lnarg <= 0)
			{
				r = -numeric_limits<double>::max();
#ifdef _DEBUG
				cerr << "Warning: Negative log in BP!\n";
#endif
//				exit(-1);	// Negative log
			}
			else
				// Update R.
				r = log(lnarg);
		}
	}
};

// Calculate the minsum decoding pi and min terms without exclusion.
struct functor_r_offms_pi {
	static double min0, min1;	// The lowest and second-lowest minima, respectively
	static int pi;					// The pi term, without exclusion.
	static inline void callback(double &q) {
		double qv = q;
		if (qv < 0)
			pi = -pi;				// This effectively does the sign function.
		qv = fabs(qv);
		if (min0 >= qv)				// Update both minima.
		{
			min1 = min0;
			min0 = qv;
		}
		else if (min1 > qv)			// Update the second-lowest minimum.
			min1 = qv;
	}
};
double functor_r_offms_pi::min0;
double functor_r_offms_pi::min1;
int functor_r_offms_pi::pi;

// Updates the R matrix with minsum decoding. Performs exclusion.
struct functor_r_offms_update {
	static inline void callback(double &r, double &q) {
		int pir = functor_r_offms_pi::pi;		// The pi term to use.

		const double qv = q;				// The Q term to use.
		// Perform exclusion on the pi term.
		if (qv < 0)
			pir = -pir;
		// Perform exclusion on the min term.
		const double qvmin = (fabs(qv) == functor_r_offms_pi::min0) ?
				functor_r_offms_pi::min1 : functor_r_offms_pi::min0;

		// Offset min sum calculation for r^(i)_(m,n)
		r = pir * max((double)0.0, qvmin - beta);
	}
};

// Calculates the sigma term without exclusion, for use in updating Q and L
// during decoding.
struct functor_sigmar {
	static double rsigma;
	static inline void callback(double &r) {
		rsigma += r;
	}
};
double functor_sigmar::rsigma;

// Update Q during decoding. Performs exclusion.
struct functor_updateq {
	// The value of Q for all elements in this column at iteration 0
	static double q0;

	static inline void callback(double &q, double &r) {
		// Performs exclusion and sets Q.
		q = q0 + functor_sigmar::rsigma - r;
	}
};
double functor_updateq::q0;


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

void setSnrDB(double snrdb)
{
	cout << "Setting the SNR to " << snrdb << " dB...\n";
	// Calculate the linear SNR and sigma from an SNR in dB
	snr = pow((double)10.0, (double)snrdb/10);
	sigma = pow((double)2.0*RATE*snr, (double)-0.5);
}


template <int valMax, int valSection>
void outputHistogram(
	Histogram<NBUCKETS, valMax, valSection, NBLOCKS> (&hists)[NSNRS][IMAX],
	const char *const name)
{
	// Output the error histograms.
	ofstream fhist;
	fhist << setprecision(10);

	// Error surface histogram
	// x - error buckets
	// y - iterations
	// z - frequency
	string filename = "hist_err_";
	filename += name;
	filename += ".tsv";
	fhist.open(filename.c_str());
	for (int i = 0; i < IMAX; i++)
	{
		for (int b = 0; b < NBUCKETS; b++)
			fhist << hists[DEFAULTSNR][i].getNormalizedFreq(b) << '\t';
		fhist << '\n';
	}
	fhist.close();

	// SNR surface histogram at zero error
	// x - SNR
	// y - iterations
	// z - frequency
	filename = "hist_snr_";
	filename += name;
	filename += ".tsv";
	fhist.open(filename.c_str());
	for (int i = 0; i < IMAX; i++)
	{
		for (snrindex = 0; snrindex < NSNRS; snrindex++)
			fhist << hists[snrindex][i].getNormalizedFreq(0) << '\t';
		fhist << '\n';
	}
	fhist.close();

	// SNR surface histogram at maximum iteration
	// x - error buckets
	// y - SNR
	// z - frequency
	filename = "hist_maxiter_";
	filename += name;
	filename += ".tsv";
	fhist.open(filename.c_str());
	for (int b = 0; b < NBUCKETS; b++)
	{
		for (snrindex = 0; snrindex < NSNRS; snrindex++)
			fhist << hists[snrindex][IMAX-1].getNormalizedFreq(b) << '\t';
		fhist << '\n';
	}
	fhist.close();
	
	// Giant 4D slice histogram
	// x - iterations
	// y - SNR
	// z - error buckets
	// size, colour - frequency
	filename = "hist_slice_";
	filename += name;
	filename += ".tsv";
	fhist.open(filename.c_str());
	for (snrindex = 0; snrindex < NSNRS; snrindex++)
	{
		for (int i = 0; i < IMAX; i++)
		{
			for (int b = 0; b < NBUCKETS; b++)
				fhist << hists[snrindex][i].getNormalizedFreq(b) << '\t';
			fhist << '\n';
		}
	}
	fhist.close();
}



void execute()
{

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

	// The decode method loop
	for (method = (DecodeMethod)0; method < ndecodes; ((int&)method)++)
	{
		cout << "Decoding using " << decodeNames[method] << " method...\n";

		// The SNR loop
		for (snrindex = 0; snrindex < NSNRS; snrindex++)
		{
			// Set the signal-to-noise ratio
			setSnrDB(snrs[snrindex]);

			int nerrs = 0;	// The number of block errors

			// The block loop
			for (int b = 1; b <= NBLOCKS; b++)
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
					cout << "Block " << b << '/' << NBLOCKS
						<< "\tBLER=" << nerrs*100.0/NBLOCKS << "%    \r";
			}
			cout << '\n';
		}
		cout << '\n';
	}

	// Output the error histograms.
	cout << "Generating histogram files...\n";

	ofstream snraxis("axis_snr.tsv");
	for (snrindex = 0; snrindex < NSNRS; snrindex++)
		snraxis << snrs[snrindex] << '\n';
	snraxis.close();

	ofstream orthaxis("axis_err_orth.tsv");
	ofstream messaxis("axis_err_mess.tsv");
	for (int b = 0; b < NBUCKETS; b++)
	{
		orthaxis << (**orthhist)->getNormalizedValFloor(b) << '\n';
		messaxis << (**messhist)->getNormalizedValFloor(b) << '\n';
	}
	orthaxis.close();
	messaxis.close();

	for (method = (DecodeMethod)0; method < ndecodes; ((int&)method)++)
	{
		string name = decodeNames[method];
		name += "_orth";
		outputHistogram(orthhist[method], name.c_str());

		name = decodeNames[method];
		name += "_mess";
		outputHistogram(messhist[method], name.c_str());
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
		orthhist[method][snrindex][i].report(functor_multhxhat::nerrs);

		int diff = 0;	 // The number of x==xhat errors
		for (int j = 0; j < Z*K; j++)
			diff += mxhat[j] != ms[j];
		messhist[method][snrindex][i].report(diff);

		if (!functor_multhxhat::nerrs)
		{
			if (diff)
			{
#ifdef _DEBUG
				cerr << "Warning: False positive; " << diff << " errors.\n";
#endif
				return false;
			}
			for (++i; i < IMAX; i++)
			{
				orthhist[method][snrindex][i].report(0);
				messhist[method][snrindex][i].report(0);
			}
			return true;
		}

		if (++i > IMAX)
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
		functor_r_offms_pi::min0 = numeric_limits<double>::max();
		functor_r_offms_pi::min1 = numeric_limits<double>::max();

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
