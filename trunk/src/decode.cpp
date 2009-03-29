/*
* $URL$
* $Date$
* $Rev$
*/

#include <cmath>
#include <limits>
#include <iostream>

#include "encode.hpp"

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

const Preaching<M,N,RHO_H_Y, RHO_H_X>	H(Ha, 0);	// Half-rate Preaching matrix H
const Preaching<M,K,RHO_HS_Y,RHO_HS_X>	Hs(Ha, 0);	// Half-rate Preaching matrix H (first half)
const Preaching<M,M,RHO_HP_Y,RHO_HP_X>	Hp(Ha, K);	// Half-rate Preaching matrix H (second half, for parity)

PreachingBased<double, M,N,RHO_H_Y,RHO_H_X> mr(H);	// R matrix
PreachingBased<double, M,N,RHO_H_Y,RHO_H_X> mq(H);	// Q matrix

// Decoding matrices
double ml[N*Z];		// L column
double ml0[N*Z];	// L column (iteration 0)
bool mxhat[N*Z];	// xhat column

DecodeMethod method;

const char *const decodeNames[ndecodes] =
{
	"bp",
	"offms"
};

// Beta (for minsum decoding)
const double beta = 0.15;

///////////////////////////////////////////////////////////////////////////////
// Functors ///////////////////////////////////////////////////////////////////

// The matrix multiplication functor for decoding orthagonality verification.
struct functor_multhxhat
{
	static int nerrs;		// The number of orthagonality errors.
	static inline void callbackProduct(int y, bool p)
	{
		// There is an error every time there is a 1 in the product.
		nerrs += p;
	}
};
int functor_multhxhat::nerrs;

// The functor to set the initial values for Q.
struct functor_setq
{
	static double l;	// The L-matrix value to take.
	static inline void callback(double &q)
	{
		q = l;
	}
};
double functor_setq::l;

// Calculate the BP decoding pi term without exclusion. Keeps a cache of tanh
// values as these take a while to get.
struct functor_r_bp_pi
{
	static double pi;					// The pi term (without exclusion)
	static double tanh_cache[RHO_H_Y];	// The cache of calculated tanh values
	static int xi;						// The current index in the row

	static inline void callback(double &q)
	{
		const double tanhval = tanh(q/2.0);	// Calculate the tanh term.
		pi *= tanhval;						// Multiply it into pi.
		tanh_cache[xi++] = tanhval;			// Cache it.
	}
};
double functor_r_bp_pi::pi;
double functor_r_bp_pi::tanh_cache[RHO_H_Y];
int functor_r_bp_pi::xi;

// The functor to update the R matrix, using BP decoding. Performs exclusion.
struct functor_r_bp_update
{
	static inline void callback(double &r)
	{
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
		}
		if (pir == 1)
		{
			r = numeric_limits<double>::max();
#ifdef _DEBUG
			cerr << "Warning: Divide by 0 in BP for r!\n";
#endif
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
			}
			else
				// Update R.
				r = log(lnarg);
		}
	}
};

// Calculate the minsum decoding pi and min terms without exclusion.
struct functor_r_offms_pi
{
	static double min0, min1;	// The lowest and second-lowest minima, respectively
	static int pi;				// The pi term, without exclusion.
	static inline void callback(double &q)
	{
		double qv = q;
		if (qv < 0)
			pi = -pi;			// This effectively does the sign function.
		qv = fabs(qv);
		if (min0 >= qv)			// Update both minima.
		{
			min1 = min0;
			min0 = qv;
		}
		else if (min1 > qv)		// Update the second-lowest minimum.
			min1 = qv;
	}
};
double functor_r_offms_pi::min0;
double functor_r_offms_pi::min1;
int functor_r_offms_pi::pi;

// Updates the R matrix with minsum decoding. Performs exclusion.
struct functor_r_offms_update
{
	static inline void callback(double &r, double &q)
	{
		int pir = functor_r_offms_pi::pi;	// The pi term to use.

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
struct functor_sigmar
{
	static double rsigma;
	static inline void callback(double &r)
	{
		rsigma += r;
	}
};
double functor_sigmar::rsigma;

// Update Q during decoding. Performs exclusion.
struct functor_updateq
{
	// The value of Q for all elements in this column at iteration 0
	static double q0;

	static inline void callback(double &q, double &r)
	{
		// Performs exclusion and sets Q.
		q = q0 + functor_sigmar::rsigma - r;
	}
};
double functor_updateq::q0;

///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

void operator++(DecodeMethod &incmethod, int)
{
	incmethod = (DecodeMethod)((int)incmethod + 1);
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
		messhist[i].report(diff);

		if (!functor_multhxhat::nerrs)
		{
			if (diff)
			{
#ifdef _DEBUG
				cerr << "Warning: False positive; " << diff << " errors.\n";
#endif
				return false;
			}
			for (i++; i < IMAX; i++)
			{
				orthhist[i].report(0);
				messhist[i].report(0);
			}
			return true;
		}

		if (++i >= IMAX)
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
		mr.iterY<functor_r_bp_update>(m);
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
