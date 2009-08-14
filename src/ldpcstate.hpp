/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include "mtrand/MTRand_gaussian.hpp"
#include "preachingbased.hpp"
#include "histogram.hpp"

// Presets for half-rate
#define M 12		// Height of the unexpanded Preaching matrix
#define N 24		// Width of the unexpanded Preaching matrix

#define OUTPUT_DEBUGFILE 0	// Enable to output data to a debug file

#define NERRBUCKETS  25		// The number of error histogram buckets
#define NPERFBUCKETS 100	// The number of performance histogram buckets
#define MAXPERFTIME  2000	// The maximum performance histogram duration, in us

#define NERRS		300
#define NSNRS		10	// The number of SNRs to try
#define DEFAULTSNR	5	// The index of the default SNR.
#define IMAX		500	// The maximum number of decode iterations

// Matrix sparsity parameters
#define RHO_H_Y  7
#define RHO_H_X  6
#define RHO_HS_Y 5
#define RHO_HS_X 6
#define RHO_HP_Y 3
#define RHO_HP_X 3


// The decode method.
struct DecodeMethod
{
	enum Enum
	{
		firstMethod = 0,	// (The first method available)
		ms = 0,				// Min sum
		ms_sc,
		offms,				// Offset min sum
		offms_sc,
		nms,				// Normalized min sum
		nms_sc,
		bp,					// Belief propagation
		ndecodes			// (The number of decoding algorithms)
	};
};

inline void operator++(DecodeMethod::Enum &incmethod, int)
{
	incmethod = (DecodeMethod::Enum)((int)incmethod + 1);
}

// Data for unexpanded half-rate Preaching matrix H
extern const int Ha[M][N];

extern const Preaching<M,N,RHO_H_Y, RHO_H_X>  H;	// Unexpanded half-rate Preaching matrix H
extern const Preaching<M,K,RHO_HS_Y,RHO_HS_X> Hs;	// Half-rate Preaching matrix H (first half)
extern const Preaching<M,M,RHO_HP_Y,RHO_HP_X> Hp;	// Half-rate Preaching matrix H (second half, for parity)

extern const char *const decodeNames[DecodeMethod::ndecodes];

// The SNRs to try.
extern const double snrs[NSNRS];

// The orthagonality error and message error histograms.
// The template parameters are the number of histogram buckets, the full size
// of the data range, and the desired portion of the data range to examine.
typedef ValNormalizedHistogram<NERRBUCKETS, M*Z, (int)(M*Z*0.33)>	OrthHistType;
typedef ValNormalizedHistogram<NERRBUCKETS, N*Z, (int)(N*Z*0.06)>	MessHistType;
typedef Histogram<NPERFBUCKETS, MAXPERFTIME>						PerfHistType;

class LDPCstate
{
public:
	PreachingBased<double, M,N,RHO_H_Y,RHO_H_X> mr;	// R matrix
	PreachingBased<double, M,N,RHO_H_Y,RHO_H_X> mq;	// Q matrix

	DecodeMethod::Enum method;

	int snrindex;	// The current SNR index
	double snr;	// Signal-to-noise ratio
	double snrdb;
	double sigma;

	bool mx[N*Z];		// (col) Combination of ms and mp
	bool mxhat[N*Z];	// xhat column
	double my[N*Z];		// (col) Encoder output after AWGN

	// Set the aliases into mx
	bool (&ms)[K*Z];	// (col) Message
	bool (&mp)[M*Z];	// (col) Generated parity

	// The Gaussian distribution random number generator
	MTRand_gaussian grand;	//((unsigned long)time(0));
	// Discrete value random number generator
	MTRand_int32 irand;		//((unsigned long)~time(0));

	bool msprod[M*Z];	// Encoding verification column

///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

	// Constructor
	LDPCstate();
	
	void init(DecodeMethod::Enum methodInit,
		int snrindexInit);

	// Computer the output of the encoder
	void encode();

	// Set the parity matrix based on the message
	void setParity();

	// Compute the output of the decoder
	// Returns true if no error
	bool decode(OrthHistType (&orthhist)[IMAX],
				MessHistType (&messhist)[IMAX],
				PerfHistType &perfhist);

	// Set the initial decoder state
	void decode_initial();

	// Update the R matrix (belief propagation method)
	void rupdate_bp();

	// Update the R matrix (minsum method)
	template <DecodeMethod::Enum msMethod, bool sc>
	void rupdate_ms();

	// Update the Q and L matrices
	template <bool sc>
	void qlupdate();

	static void calculateRho();
};
