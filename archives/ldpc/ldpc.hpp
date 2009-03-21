#ifndef INCLUDED_LDPC_HPP
#define INCLUDED_LDPC_HPP

#include "automatrix.hpp"
#include "preaching.hpp"
#include "mtrand/MTRand_gaussian.hpp"

#include <cmath>

// Presets for half-rate
#define M 12	// Height of the unexpanded Preaching matrix
#define N 24	// Width of the unexpanded Preaching matrix

class LDPC
{
private:
	// Data for unexpanded half-rate Preaching matrix H
	static const int Ha[M][N];

	static const Preaching<M, N, 7> H;  // Unexpanded half-rate Preaching matrix H
	static const Preaching<M, K, 6> Hs; // Unexpanded half-rate Preaching matrix H (first half)
	static const Preaching<M, M, 3> Hp; // Unexpanded half-rate Preaching matrix H (second half, for parity)

	// Encoding matrices
	static bool (&ms)[K*Z];						// (col) Message
	static bool (&mp)[M*Z];						// (col) Generated parity
	static Automatrix1<bool, N*Z> mx;			// (col) Combination of ms and mp

	// Encoding verification column
	static Automatrix1<bool, M*Z> msprod;

	// Encoded message after noise channel
	static Automatrix1<long double, N*Z> my;	// (col) Encoder output after AWGN

	// Decoding matrices
	static Automatrix2<long double, M*Z, N*Z> mr;
	static Automatrix2<long double, N*Z, M*Z> mq;
	static Automatrix2<long double, N*Z, M*Z> mq0;
	static Automatrix1<long double, N*Z> ml;
	static Automatrix1<long double, N*Z> ml0;
	static Automatrix1<bool, N*Z> mxhat;


	// The Gaussian distribution random number generator
	MTRand_gaussian grand;	
	// Discrete value random number generator
	MTRand_int32 irand;


	// Constants for AWGN calculation
	static const double R;
	double snr;
	double snrdb;
	double sigma;

	// Constants for decoding
	unsigned imax;


public:
	LDPC();

	// Run the simulation
	void execute();

	// Computer the output of the encoder
	void encode();

	// Compute the output of the decoder
	// Returns 0 if no error, or 1 if there was an error
	unsigned decode();
	
	void rupdate_bp();
	
	void rupdate_offms();






	// Functor to find the sum over elements in ms for a row,
	// iterating in x
	struct functor_summsy
	{
		static bool sum;
		static inline void callback(unsigned y, unsigned x)
		{
			sum ^= ms[x];
		}
	};

	struct functor_multhpp
	{
		static bool success;
		static inline void callback(unsigned y, bool p)
		{
			success &= msprod[y] == p;
		}
	};

	struct functor_multhxhat
	{
		static unsigned nerrs;
		static inline void callback(unsigned y, bool p)
		{
			nerrs += p;
		}
	};

	struct functor_r_bp
	{
		static long double pi;
		static inline void callback(unsigned y, unsigned x)
		{
			pi *= tanh(mq[x][y]/2); // WRONG?
		}
	};
	
	struct functor_r_offms
	{
		static int pi;
		static long double min0, min1;
		static inline void callback(unsigned y, unsigned x)
		{
			long double qv = mq[y][x];
			pi *= (qv > 0) ? 1 : -1;
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

	struct functor_sigmar
	{
		static long double sigma;
		static inline void callback(unsigned y, unsigned x)
		{
			sigma += mr[y][x];
		}
	};




private:
	
	// Set the parity matrix based on the message
	void setParity();
};

#endif
