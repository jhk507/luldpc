#pragma once

#include "automatrix.hpp"
#include "preachingbased.hpp"
#include "mtrand/MTRand_gaussian.hpp"

#include <cmath>
#include <fstream>

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
	static bool (&ms)[K*Z];				// (col) Message
	static bool (&mp)[M*Z];				// (col) Generated parity
	static Automatrix1<bool, N*Z> mx;	// (col) Combination of ms and mp

	// Encoding verification column
	static Automatrix1<bool, M*Z> msprod;

	// Encoded message after noise channel
	static Automatrix1<long double, N*Z> my;	// (col) Encoder output after AWGN

	// Decoding matrices
	static PreachingBasedR<long double, M, N> mr; //check node to variable node
	static PreachingBasedQ<long double, M, N> mq; //variable node to check node
	static PreachingBasedQ<long double, M, N> mq0; //variable node to check node at 0th iteration
	static Automatrix1<long double, N*Z> ml; // LLR update
	static Automatrix1<long double, N*Z> ml0; // ln^(0)
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

	static const double beta;

	std::ofstream debugfile;


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


	struct functor_setq
	{
		static long double l;
		static inline void callbackX(unsigned X, unsigned i)
		{
			//mq0.Vc[y][i] = l;
			//mq.Vc[y][i] = l;
		}
	};

	struct functor_updateq
	{
		static inline void callbackY()
		{
		}
	};


	// Functor to find the sum over elements in ms for a row,
	// iterating in x
	struct functor_summsy
	{
		static bool sum;
		static inline void callback(unsigned y, unsigned x, unsigned i)
		{
			sum ^= ms[x]; //XOR
		}
	};

	struct functor_multhpp
	{
		static bool success;
		static inline void callback(unsigned y, bool p)
		{
			success &= msprod[y] == p;  //checks if message sum product = parity SP
		}
	};

	struct functor_multhxhat
	{
		static unsigned nerrs;
		static inline void callback(unsigned y, bool p)
		{
			nerrs += p; //error
		}
	};

	struct functor_r_bp
	{
		static long double pi;
		static inline void callback(unsigned y, unsigned x, unsigned i)
		{
			// pi *= tanh(mq[x][y]/2);
			pi *= tanh(mq.Vc[y][i]/2);
		}
	};


	struct functor_r_offms
	{
		static int pi;
		static long double min0, min1;
		static inline void callback(unsigned y, unsigned x, unsigned i)
		{
			// long double qv = mq[y][x];
			long double qv = mq.Vc[y][i];
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
		static inline void callback(unsigned y, unsigned x, unsigned i)
		{
			sigma += mr.Vc[x][i]; //update
		}
	};

private:

	// Set the parity matrix based on the message
	void setParity();

	template <int X, typename MatrixType>
	static void outputLargeMatrix1(MatrixType &matrix)
	{
		for (int x = 0; x < X*Z; x++)
		{
			debugfile << matrix[x] << '\t';
			if (x%Z == Z-1)
				debugfile << '\t';
		}
		debugfile << endl;
	}

	template <int X>
	static void outputLargeMatrix1<X, bool(&)[X*Z]>(bool (&matrix)[X*Z])
	{
		bool *pm = matrix;
		for (int x = 0; x < X; x++)
		{
			for (int z = 0; z < Z; z++, pm++)
				debugfile << *pm;
			debugfile << ' ';
		}
		debugfile << endl;
	}

	template <int Y, int X, typename MatrixType>
	static void outputLargeMatrix2(MatrixType &matrix)
	{
		for (int y = 0; y < Y; y++)
		{
			for (int x = 0; x < X; x++)
			{
				debugfile << matrix[y][x] << '\t';
				if (x%Z == Z-1)
					debugfile << '\t';
			}
			debugfile << endl;
			if (y%Z == Z-1)
				debugfile << endl;
		}
		debugfile << endl;
	}
};
