/*
* $URL$
* $Date$
* $Rev$
*/

#include <limits>
#include <iostream>

#include "itime.hpp"
#include "mathfun.hpp"
#include "ldpcstate.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Globals ////////////////////////////////////////////////////////////////////

// Data for unexpanded half-rate Preaching matrix H
const int Ha[M][N] =
{
#if (M == 12 && N == 24)									// rate = 1/2
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
#elif(M == 6 && N == 24)
#if (VARIANT == 'A')										// rate = 3/4, A
	{ 6,38, 3,93,-1,-1,-1,30,70,-1,86,-1,37,38, 4,11,-1,46,48, 0,-1,-1,-1,-1},
	{62,94,19,84,-1,92,78,-1,15,-1,-1,92,-1,45,24,32,30,-1,-1, 0, 0,-1,-1,-1},
	{71,-1,55,-1,12,66,45,79,-1,78,-1,-1,10,-1,22,55,70,82,-1,-1, 0, 0,-1,-1},
	{38,61,-1,66, 9,73,47,64,-1,39,61,43,-1,-1,-1,-1,95,32, 0,-1,-1, 0, 0,-1},
	{-1,-1,-1,-1,32,52,55,80,95,22, 6,51,24,90,44,20,-1,-1,-1,-1,-1,-1, 0, 0},
	{-1,63,31,88,20,-1,-1,-1, 6,40,56,16,71,53,-1,-1,27,26,48,-1,-1,-1,-1, 0}
#elif (VARIANT == 'B')										// rate = 3/4, B
	{-1,81,-1,28,-1,-1,14,25,17,-1,-1,85,29,52,78,95,22,92, 0, 0,-1,-1,-1,-1},
	{42,-1,14,68,32,-1,-1,-1,-1,70,43,11,36,40,33,57,38,24,-1, 0, 0,-1,-1,-1},
	{-1,-1,20,-1,-1,63,39,-1,70,67,-1,38, 4,72,47,29,60, 5,80,-1, 0, 0,-1,-1},
	{64, 2,-1,-1,63,-1,-1, 3,51,-1,81,15,94, 9,85,36,14,19,-1,-1,-1, 0, 0,-1},
	{-1,53,60,80,-1,26,75,-1,-1,-1,-1,86,77, 1, 3,72,60,25,-1,-1,-1,-1, 0, 0},
	{77,-1,-1,-1,15,28,-1,35,-1,72,30,68,85,84,26,64,11,89, 0,-1,-1,-1,-1, 0}
#endif
#elif (M == 8 && N == 24)
#if (VARIANT == 'A')										// rate = 2/3, A
	{ 3, 0,-1,-1, 2, 0,-1, 3, 7,-1, 1, 1,-1,-1,-1,-1, 1, 0,-1,-1,-1,-1,-1,-1},
	{-1,-1, 1,-1,36,-1,-1,34,10,-1,-1,18, 2,-1, 3, 0,-1, 0, 0,-1,-1,-1,-1,-1},
	{-1,-1,12, 2,-1,15,-1,40,-1, 3,-1,15,-1, 2,13,-1,-1,-1, 0, 0,-1,-1,-1,-1},
	{-1,-1,19,24,-1, 3, 0,-1, 6,-1,17,-1,-1,-1, 8,39,-1,-1,-1, 0, 0,-1,-1,-1},
	{20,-1, 6,-1,-1,10,29,-1,-1,28,-1,14,-1,38,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1},
	{-1,-1,10,-1,28,20,-1,-1, 8,-1,36,-1, 9,-1,21,45,-1,-1,-1,-1,-1, 0, 0,-1},
	{35,25,-1,37,-1,21,-1,-1, 5,-1,-1, 0,-1, 4,20,-1,-1,-1,-1,-1,-1,-1, 0, 0},
	{-1, 6, 6,-1,-1,-1, 4,-1,14,30,-1, 3,36,-1,14,-1, 1,-1,-1,-1,-1,-1,-1, 0}
#elif (VARIANT == 'B')										// rate = 2/3, B
	{ 2,-1,19,-1,47,-1,48,-1,36,-1,82,-1,47,-1,15,-1,95, 0,-1,-1,-1,-1,-1,-1},
	{-1,69,-1,88,-1,33,-1, 3,-1,16,-1,37,-1,40,-1,48,-1, 0, 0,-1,-1,-1,-1,-1},
	{10,-1,86,-1,62,-1,28,-1,85,-1,16,-1,34,-1,73,-1,-1,-1, 0, 0,-1,-1,-1,-1},
	{-1,28,-1,32,-1,81,-1,27,-1,88,-1, 5,-1,56,-1,37,-1,-1,-1, 0, 0,-1,-1,-1},
	{23,-1,29,-1,15,-1,30,-1,66,-1,24,-1,50,-1,62,-1,-1,-1,-1,-1, 0, 0,-1,-1},
	{-1,30,-1,65,-1,54,-1,14,-1, 0,-1,30,-1,74,-1, 0,-1,-1,-1,-1,-1, 0, 0,-1},
	{32,-1, 0,-1,15,-1,56,-1,85,-1, 5,-1, 6,-1,52,-1, 0,-1,-1,-1,-1,-1, 0, 0},
	{-1, 0,-1,47,-1,13,-1,61,-1,84,-1,55,-1,78,-1,41,95,-1,-1,-1,-1,-1,-1, 0}
#endif
#elif (M == 4 && N == 24)									// rate = 5/6
	{ 1,25,55,-1,47, 4,-1,91,84, 8,86,52,82,33, 5, 0,36,20, 4,77,80, 0,-1,-1},
	{-1, 6,-1,36,40,47,12,79,47,-1,41,21,12,71,14,72, 0,44,49, 0, 0, 0, 0,-1},
	{51,81,83, 4,67,-1,21,-1,31,24,91,61,81, 9,86,78,60,88,67,15,-1,-1, 0, 0},
	{50,-1,50,15,-1,36,13,10,11,20,53,90,29,92,57,30,84,92,11,66,80,-1,-1, 0}
#endif
};

const Preaching<M,N,RHO_H_Y, RHO_H_X>	H(Ha, 0);	// Half-rate Preaching matrix H
const Preaching<M,K,RHO_HS_Y,RHO_HS_X>	Hs(Ha, 0);	// Half-rate Preaching matrix H (first half)
const Preaching<M,M,RHO_HP_Y,RHO_HP_X>	Hp(Ha, K);	// Half-rate Preaching matrix H (second half, for parity)

const char *const decodeNames[DecodeMethod::ndecodes] =
{
	"ms",
	"ms_sc",
	"offms",
	"offms_sc",
	"nms",
	"nms_sc",
	"v_off_ms",
	"bp"
};

// Alpha (for normalized minsum decoding)
const double alpha   = 0.8;
const double alphasc = 0.92;
// Beta (for minsum decoding)
const double beta   = 0.15;
const double betasc = 0.08;
const double betav  = 0.15; 

///////////////////////////////////////////////////////////////////////////////
// Functors ///////////////////////////////////////////////////////////////////

// The matrix multiplication functor for decoding orthagonality verification.
struct functor_multhxhat
{
	int nerrs;		// The number of orthagonality errors.
	inline void callbackProduct(int y, bool p)
	{
		// There is an error every time there is a 1 in the product.
		nerrs += p;
	}
};

// The functor to set the initial values for Q.
struct functor_setq
{
	double l;	// The L-matrix value to take.
	inline void callback(double &q)
	{
		q = l;
	}
};

// Calculate the BP decoding pi term without exclusion. Keeps a cache of tanh
// values as these take a while to get.
struct functor_r_bp_pi
{
	double pi;					// The pi term (without exclusion)
	double tanh_cache[RHO_H_Y];	// The cache of calculated tanh values
	double *pcache;				// The current index in the row

	inline void callback(double &q)
	{
		const double tanhval = tanh/*app*/(q/2.0);	// Calculate the tanh term.
		pi *= tanhval;							// Multiply it into pi.
		*pcache = tanhval;						// Cache it.
		pcache++;
	}
};

// The functor to update the R matrix, using BP decoding. Performs exclusion.
struct functor_r_bp_update
{
	double pir0;			// The starting pi term to use.
	double *tanh_cache;		// The appropriate cached tanh to start.

	inline void callback(double &r)
	{
		double pir = pir0;
		
		const double tanhr = *tanh_cache;
		tanh_cache++;

		if (!tanhr)
		{
#ifdef _DEBUG
			cerr << "Warning: Divide by 0; tanhr=0 in BP!\n";
#endif
			// lim p->inf ln((1+p)/(1-p))
			//          = ln(-1)
			r = -1e20; // -numeric_limits<double>::max();
			return;
		}

		// Perform exclusion.
		pir /= tanhr;

		if (pir == 1)
		{
#ifdef _DEBUG
			cerr << "Warning: Divide by 0; pir=1 in BP!\n";
#endif
			// lim p->1 ln((1+p)/(1-p))
			//        = ln(2/0)
			//        = inf
			r = 1e20; // numeric_limits<double>::max();
			return;
		}

		// Calculate the argument for ln.
		const double lnarg = (1+pir)/(1-pir);

		if (lnarg <= 0)
		{
#ifdef _DEBUG
			cerr << "Warning: Negative log in BP!\n";
#endif
			r = -1e20; // -numeric_limits<double>::max();
			return;
		}

		// Update R.
		r = log(lnarg);
	}
};

// Calculate the minsum decoding pi and min terms without exclusion.
struct functor_r_ms_pi
{
	double min0, min1;	// The lowest and second-lowest minima, respectively
	int pi;				// The pi term, without exclusion.
	inline void callback(double &q)
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

// Updates the R matrix with minsum decoding. Performs exclusion.
template <DecodeMethod::Enum msMethod, bool sc>
struct functor_r_ms_update
{
	int pir0;
	double min0, min1;
	inline void callback(double &r, double &q)
	{
		int pir = pir0;	// The pi term to use.

		const double qv = q;				// The Q term to use.
		// Perform exclusion on the pi term.
		if (qv < 0)
			pir = -pir;
		// Perform exclusion on the min term.
		const double qvmin = (fabs(qv) == min0) ? min1 : min0;

		// Offset min sum calculation for r^(i)_(m,n)
		r = pir;
		if (msMethod == DecodeMethod::offms)
			r *= max((double)0.0, qvmin - (sc ? betasc : beta));
		else
			r *= qvmin;

		if (msMethod == DecodeMethod::nms)
			r *= sc ? alphasc : alpha;
	}
};

// Calculates the sigma term without exclusion, for use in updating Q and L
// during decoding.
struct functor_sigmar
{
	double rsigma;
	inline void callback(double &r)
	{
		rsigma += r;
	}
};

// Update Q during decoding. Performs exclusion.
template <DecodeMethod::Enum method, bool sc>
struct functor_updateq
{
	// The value of Q for all elements in this column at iteration 0
	double q0;
	double rsigma;

	inline void callback(double &q, double &r)
	{
		// Performs exclusion and sets Q.
		const double qnew = q0 + rsigma - r;
		if (method == DecodeMethod::v_off_ms)
		{
			const bool qsgn = qnew < 0;
			double qmag = fabs(qnew);
			if (qmag > betav)
			{
				qmag -= betav;
				if (qsgn) qmag = -qmag;
				q = qmag;
			}
			else
			{
				if (r != 0 && (r < 0) != qsgn)
					q = 0;
			}
		}
		else
		{
			if (sc)
			{
				if (q != 0 && ((q >= 0) != (qnew >= 0)))
					q = 0;
				else
					q = qnew;
			}
			else
				q = qnew;
		}
	}
};

///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

// Do block-by-block profiling.
// Profiler profs[8];

bool LDPCstate::decode(
	OrthHistType (&orthhist)[IMAX],
	MessHistType (&messhist)[IMAX],
	PerfHistType &perfhist)
{
	// Set the initial state of the decoder
	decode_initial();

	functor_multhxhat funcverify;

	// Iterative decoding
	for (iter = 0; ; )
	{
#if OUTPUT_DEBUGFILE
		debugfile << "Before iteration " << iter << ":" << endl;

		debugfile << "l:" << endl;
		outputLarge<N,Z>(ml, debugfile);

		debugfile << "q:" << endl;
		mq.output(debugfile);

		debugfile.flush();

		if (iter >= 1)
			return true;
#endif

		ITime perf;

		// Update the R matrix
		switch (method)
		{
		case DecodeMethod::ms:			rupdate_ms<DecodeMethod::ms,	false>();	break;
		case DecodeMethod::ms_sc:		rupdate_ms<DecodeMethod::ms, 	true >();	break;
		case DecodeMethod::offms:		rupdate_ms<DecodeMethod::offms,	false>();	break;
		case DecodeMethod::offms_sc:	rupdate_ms<DecodeMethod::offms,	true >();	break;
		case DecodeMethod::nms:			rupdate_ms<DecodeMethod::nms,	false>();	break;
		case DecodeMethod::nms_sc:		rupdate_ms<DecodeMethod::nms,	true >();	break;
		case DecodeMethod::v_off_ms:	rupdate_ms<DecodeMethod::v_off_ms,	false>();	break; 
		case DecodeMethod::bp:			rupdate_bp();								break;
		}

#if OUTPUT_DEBUGFILE
		debugfile << "After iteration " << iter << ":" << endl;

		debugfile << "r:" << endl;
		mr.output(debugfile);

		debugfile.flush();
#endif

		// Update the Q and L matrices
		switch (method)
		{
		case DecodeMethod::ms:
		case DecodeMethod::offms:
		case DecodeMethod::nms:
			qlupdate<DecodeMethod::ms, false>();
			break;
		case DecodeMethod::ms_sc:
		case DecodeMethod::offms_sc:
		case DecodeMethod::nms_sc:
			qlupdate<DecodeMethod::ms, true>();
			break;
		case DecodeMethod::v_off_ms:
			qlupdate<DecodeMethod::v_off_ms, false>();
			break;
		case DecodeMethod::bp:
			qlupdate<DecodeMethod::bp, false>();
		}

		funcverify.nerrs = 0;
		H.multColCallback(funcverify, mxhat);
		orthhist[iter].report(funcverify.nerrs);

		int diff = 0;	 // The number of x==xhat errors
		for (int j = 0; j < Z*K; j++)
			diff += mxhat[j] != ms[j];
		messhist[iter].report(diff);

		perfhist.report((int)(perf.get()*1e6));

		if (!funcverify.nerrs)
		{
			if (diff)
			{
#ifdef _DEBUG
				cerr << "Warning: False positive; " << diff << " errors.\n";
#endif
				return false;
			}
			for (int i = ++iter; i < IMAX; i++)
			{
				orthhist[i].report(0);
				messhist[i].report(0);
			}
			return true;
		}

		if (++iter >= IMAX)
			return false;
	}
}

void LDPCstate::decode_initial()
{
	functor_setq funcq;

	// Set initial state
	for (int n = 0; n < N*Z; n++)
	{
		if (method == DecodeMethod::bp)
			my[n] *= 2.0/sigma/sigma; // Required for BP algorithm

		funcq.l = my[n];
		mq.iterX(funcq, n);
	}
}

void LDPCstate::rupdate_bp()
{
	functor_r_bp_pi funcpi;
	functor_r_bp_update funcup;
	// Update mr
	for (int m = 0; m < Z*M; m++)
	{
		funcpi.pi = 1;
		funcpi.pcache = funcpi.tanh_cache;
		mq.iterY(funcpi, m);

		// Reset the row index
		funcup.tanh_cache = funcpi.tanh_cache;
		funcup.pir0 = funcpi.pi;
		mr.iterY(funcup, m);
	}
}

template <DecodeMethod::Enum msMethod, bool sc>
void LDPCstate::rupdate_ms()
{
	// Update mr
	// OFF-MS method
	functor_r_ms_pi funcpi;
	functor_r_ms_update<msMethod, sc> funcup;
	for (int m = 0; m < Z*M; m++)
	{
		funcpi.pi = 1;	// Multiplicative identity

		// Min function identities
		// Assume a very large value so that it may be overwritten on the first
		// iteration.
		funcpi.min0 = numeric_limits<double>::max();
		funcpi.min1 = numeric_limits<double>::max();
		mq.iterY(funcpi, m);

		funcup.min0 = funcpi.min0;
		funcup.min1 = funcpi.min1;
		funcup.pir0 = funcpi.pi;
		mr.iterY2(funcup, m,mq);
	}
}

template <DecodeMethod::Enum msMethod, bool sc>
void LDPCstate::qlupdate()
{
	functor_sigmar funcr;
	functor_updateq<msMethod, sc> funcq;
	// Update mq and ml
	for (int n = 0; n < Z*N; n++)
	{
		funcr.rsigma = 0;
		mr.iterX(funcr, n);

		funcq.q0 = my[n];
		funcq.rsigma = funcr.rsigma;
		mq.iterX2(funcq, n,mr);

		const double ml = funcq.q0 + funcr.rsigma;
		mxhat[n] = ml < 0; // Hard decision
	}
}
