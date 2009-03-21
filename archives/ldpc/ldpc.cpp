#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>

#include "ldpc.hpp"

using namespace std;

// Data for unexpanded half-rate Preaching matrix H
const int LDPC::Ha[M][N] =
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

const Preaching<M,N,7> LDPC::H(Ha, 0);
const Preaching<M,K,6> LDPC::Hs(Ha, 0); // Unexpanded half-rate Preaching matrix H (first half)
const Preaching<M,M,3> LDPC::Hp(Ha, K); // Unexpanded half-rate Preaching matrix H (second half, for parity)

Automatrix1<bool, N*Z> LDPC::mx;
Automatrix1<long double, N*Z> LDPC::my;

// Set the aliases into mx
bool (&LDPC::ms)[K*Z] = (bool(&)[K*Z])*mx.getData(0);
bool (&LDPC::mp)[M*Z] = (bool(&)[M*Z])*mx.getData(K*Z);

Automatrix1<bool, M*Z> LDPC::msprod;

Automatrix2<long double, M*Z, N*Z> LDPC::mr;
Automatrix2<long double, N*Z, M*Z> LDPC::mq;
Automatrix2<long double, N*Z, M*Z> LDPC::mq0;
Automatrix1<long double, N*Z> LDPC::ml;
Automatrix1<long double, N*Z> LDPC::ml0;
Automatrix1<bool, N*Z> LDPC::mxhat;

const double LDPC::R = 0.5;

bool        LDPC::functor_summsy::sum;
bool        LDPC::functor_multhpp::success;
long double LDPC::functor_r_bp::pi;
long double LDPC::functor_sigmar::sigma;
unsigned    LDPC::functor_multhxhat::nerrs;
int         LDPC::functor_r_offms::pi;
long double LDPC::functor_r_offms::min0;
long double LDPC::functor_r_offms::min1;




LDPC::LDPC() :
	grand((unsigned long)time(0)),
	irand((unsigned long)~time(0))
{
	cout << "Enter signal to noise ratio (dB): ";
	cin >> snrdb;
	snr = pow(10.0,snrdb/10);
	sigma = pow(2*R*snr, -0.5);
	cout << "Enter maximum decoding iteration count: ";
	cin >> imax;
}



void LDPC::execute()
{
	unsigned nerrs = 0;	// The number of block errors

	// The main block loop
	for (unsigned b = 1; nerrs < 50; b++)
	{
		// Encode
		encode();

		// Output y
//		ofstream yout("yout.csv");
//		for (unsigned i = 0; i < N*Z; i++)
//			yout << my[i] << endl;
			

		// Decode
		nerrs += decode();

		cout << "Block " << b << ": " << nerrs << " errors, BLER=" << 100.0*nerrs/b << '%' << endl;
	}
}


void LDPC::encode()
{
	// Generate the random bits of the message
	for (unsigned m = 0; m < K*Z; m++)
		ms[m] = irand() & 1;

	// Encode
#ifdef _DEBUG
	cout << "Encoding..." << endl;
#endif

	// Get the parity bits
	setParity();

#ifdef _DEBUG
	// Double-check that the encoding succeeded
	functor_multhpp::success = true;
	Hs.multCol(ms, msprod);
	Hp.multCol<functor_multhpp>(mp);
	cout << "Encoder check " << (functor_multhpp::success ? "passed." : "failed.") << endl;
#endif

	for (unsigned n = 0; n < N*Z; n++)
		// Perform BPSK and AWGN addition
		my[n] = (mx[n] ? -1 : 1) + grand()*sigma;	// Uses ternary operator for BPSK
}

// Set the parity matrix based on the message
void LDPC::setParity()
{
	// Magic number for left side of v(0) determination
	// With our matrix, this element is ZERO therefore there is NO SHIFT NEEDED
	// const int xshift = 0;

	// Determine v(0)
	for (unsigned mi = 0; mi < Z; mi++) // Iterate over the index of v0
	{
		functor_summsy::sum = 0;
		for (unsigned m = mi; m < Z*M; m += Z)	// Iterate over m for whole H matrix
			Hs.iterY<functor_summsy>(m);		// Effectively iterate over n
		//mp[(mi+xshift)%Z] = sum;
		mp[mi] = functor_summsy::sum;
	}

/*
	// v(0) alternate method check
	Automatrix1<bool, M*Z> mpcheck;
	for (unsigned i = 0; i < Z; i++)
		mpcheck[i] = 0;
	for (unsigned m = 0; m < M; m++)
		for (unsigned n = 0; n < K; n++)
			for (unsigned mi = 0; mi < Z; mi++)
				for (unsigned ni = 0; ni < Z; ni++)
					if (Hm.at(Z*m+mi,Z*n+ni))
						mpcheck[mi] ^= ms[n*Z+ni];
	for (unsigned i = 0; i < Z; i++)
		if (mpcheck[i] != mp[i])
			cout << "FAIL" << endl;
*/

	// Determine v(1)
	for (unsigned mi = 0; mi < Z; mi++)
	{
		functor_summsy::sum = Hp.pshift(0,0,mp,mi);	// P(i,k)v(0)
		Hs.iterY<functor_summsy>(mi);				// sigma P(i,j)u(j)
		mp[Z+mi] = functor_summsy::sum;				// p(1)
	}
/*
	// v(1) alternate method check
	for (unsigned mi = 0; mi < Z; mi++)
		for (unsigned ni = 0; ni < Z; ni++)
			if (Hp.at(mi, ni))
				mpcheck[mi+Z] = mpcheck[ni];
	for (unsigned n = 0; n < M; n++)
		for (unsigned mi = 0; mi < Z; mi++)
			for (unsigned ni = 0; ni < Z; ni++)
				if (Hm.at(mi, n*Z+ni))
					mpcheck[mi+Z] ^= ms[n*Z+ni];
	for (unsigned i = Z; i < 2*Z; i++)
		if (mpcheck[i] != mp[i])
				cout << "FAIL" << endl;
*/
	// Determine v(i)
	bool *pmp = mp+Z;	// p(i) starting at i=1
	for (unsigned i = 1; i <= M-2; i++)
	{
		for (unsigned mi = 0; mi < Z; mi++, pmp++)
		{
			// v(i) + P(i,k)v(0)
			functor_summsy::sum = *pmp ^ Hp.pshift(i,0,mp,mi);

			// sigma P(i,j)u(j)
			Hs.iterY<functor_summsy>(Z*i+mi);

			pmp[Z] = functor_summsy::sum;	// p(i+1)
		}
	}
/*
	// v(i) alternate method check
	for (unsigned i = 2; i < M; i++)
	{
		for (unsigned mi = 0; mi < Z; mi++)
		{
			mpcheck[i*Z+mi] = mpcheck[(i-1)*Z+mi];
			for (unsigned ni = 0; ni < Z; ni++)
				if (Hp.at((i-1)*Z+mi, ni))
					mpcheck[i*Z+mi] ^= mpcheck[ni];
		}

		for (unsigned n = 0; n < K; n++)
			for (unsigned mi = 0; mi < Z; mi++)
				for (unsigned ni = 0; ni < Z; ni++)
					if (Hm.at((i-1)*Z+mi,n*Z+ni))
						mpcheck[i*Z+mi] ^= ms[n*Z+ni];
	}
	for (unsigned i = 2*Z; i < M*Z; i++)
		if (mpcheck[i] != mp[i])
			cout << "FAIL" << endl;
*/
}


unsigned LDPC::decode()
{
	// Set initial state
	for (unsigned n = 0; n < N*Z; n++)
	{
		const long double l = 2.0/sigma/sigma*my[n];
		ml0[n] = l;
		ml[n] = l;
		for (unsigned m = 0; m < M*Z; m++)
		{
			mq0[n][m] = l;
			mq[n][m] = l;
		}
	}

	// Iterative decoding
	for (unsigned i = 0; ; )
	{
		// Update mr
		//rupdate_bp();
		rupdate_offms();
		
		// Update mq and ml
		for (unsigned n = 0; n < Z*N; n++)
		{
			functor_sigmar::sigma = 0;

			H.iterX<functor_sigmar>(n);

			for (unsigned m = 0; m < Z*M; m++)
			{
				mq[n][m] = mq0[n][m] + functor_sigmar::sigma;
				if (H.at(m, n))
					mq[n][m] -= mr[m][n];
			}

			ml[n] = ml0[n] + functor_sigmar::sigma;

			mxhat[n] = ml[n] < 0;
		}

		// Check that the encoding succeeded
		functor_multhxhat::nerrs = 0;
		H.multCol<functor_multhxhat>(mxhat);
		
		unsigned diff = 0;
		for (int j = 0; j < Z*K; j++)
			diff += mxhat[j] != ms[j];

#ifdef _DEBUG
		cout << "Iteration " << setw(3) << i << ": "
			<< setw(4) << functor_multhxhat::nerrs << " H*xhat errors, "
			<< setw(4) << diff << " x==xhat errors."
			<< endl;
#endif

		if (!functor_multhxhat::nerrs)
			return 0;

		if (++i > imax)
		{
			cout << "Maximum iteration reached; error." << endl;
			return 1;
		}
	}
}

void LDPC::rupdate_bp()
{
	// Update mr
	for (unsigned m = 0; m < Z*M; m++)
	{
		// Do the graph iteration to calculate the pi term without exclusion
		functor_r_bp::pi = 1;

		H.iterY<functor_r_bp>(m);

		for (unsigned n = 0; n < Z*N; n++)
		{
			long double pir = functor_r_bp::pi;
			if ((1+pir)/(1-pir)<0)
			{
				n=n; // Explode!
			}
			if (H.at(m, n))
			{
				const long double tanhr = tanh(mq[n][m]/2);
				if (tanhr)
				{
					pir /= tanhr;
					if ((1+pir)/(1-pir)<0)
					{
						n=n; // Explode!
					}
				}
				else
				{
					n=n; // Explode!
				}
			}
			if (pir == 1)
				mr[m][n] = numeric_limits<long double>::max(); // Explode!
			else
				mr[m][n] = log((1+pir)/(1-pir));
		}
	}
}

void LDPC::rupdate_offms()
{
	// Update mr
	// OFF-MS method
	for (int m = 0; m < Z*M; m++)
	{
		functor_r_offms::pi = 1;	// Multiplicative identity
		// Min function identities
		functor_r_offms::min0 = numeric_limits<long double>::max();
		functor_r_offms::min1 = functor_r_offms::min0;

		// Do the graph iteration to calculate the pi and min
		// terms without exclusion
		H.iterY<functor_r_offms>(m);

		for (int n = 0; n < Z*N; n++)
		{
			long double qv = mq[n][m];
			long double qvmin;
			// Perform exclusion if necessary
			if (H.at(m, n))
			{
				functor_r_offms::pi /= (qv > 0) ? 1 : -1;
				qvmin = (fabs(qv) == functor_r_offms::min0) ? functor_r_offms::min1 : functor_r_offms::min0;
			}
			else
				qvmin = functor_r_offms::min0;

			mr[m][n] = functor_r_offms::pi * max((long double)0.0, qvmin);
		}
	}
}
