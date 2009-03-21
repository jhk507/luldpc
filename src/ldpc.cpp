#include <ctime>
#include <iostream>
#include <iomanip>
#include <limits>

#define OUTPUT_DEBUGFILE 0

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

PreachingBasedR<long double, M, N> LDPC::mr(H);
PreachingBasedQ<long double, M, N> LDPC::mq(H);
PreachingBasedQ<long double, M, N> LDPC::mq0(H);
Automatrix1<long double, N*Z> LDPC::ml;
Automatrix1<long double, N*Z> LDPC::ml0;
Automatrix1<bool, N*Z> LDPC::mxhat;

const double LDPC::R = 0.5; //rate 
const double LDPC::beta = 0.15; // Beta for minsum decoding

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
	irand((unsigned long)~time(0)),
	debugfile("debugfile.tsv")
{
	cout << "Enter signal to noise ratio (dB): ";
//	cin >> snrdb;
	snrdb = 1.5;
	snr = pow(10.0,snrdb/10);
	sigma = pow(2*R*snr, -0.5);
	cout << "Enter maximum decoding iteration count: ";
//	cin >> imax;
	imax = 50;
}

void LDPC::execute()
{
	unsigned nerrs = 0;	// The number of block errors

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
	for (int m = 0; m < M*Z; m++)
	{
		for (int n = 0; n < N*Z; n++)
		{
			debugfile << H.at(m,n) << '\t';
			if (n%Z == Z-1)
				debugfile << '\t';
		}
		debugfile << endl;
		if (m%Z == Z-1)
			debugfile << endl;
	}
	debugfile << endl;
	debugfile.flush();
#endif

	// The main block loop
	for (int b = 1; nerrs < 50; b++)
	{
		// Encode
		encode();

#if OUTPUT_DEBUGFILE
		debugfile << "Message:" << endl;
		for (int i = 0; i < K*Z; i++)
			debugfile << ms[i] << '\t';
		debugfile << endl;

		debugfile << "Encoded parity bits:" << endl;
		for (int i = 0; i < M*Z; i++)
			debugfile << mp[i] << '\t';
		debugfile << endl;
#endif

		// Decode
		nerrs += decode();

		cout << "Block " << b << ": " << nerrs << " errors, BLER=" << 100.0*nerrs/b << '%' << endl;
		debugfile << "Block " << b << ": " << nerrs << " errors, BLER=" << 100.0*nerrs/b << '%' << endl;
#if OUTPUT_DEBUGFILE
		break;
#endif
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
#if OUTPUT_DEBUGFILE
	debugfile << "Encoding..." << endl;
#endif
#endif

	// Get the parity bits
	setParity();

#ifdef _DEBUG
	// Double-check that the encoding succeeded
	functor_multhpp::success = true;
	Hs.multCol(ms, msprod);
	Hp.multCol<functor_multhpp>(mp);
	cout << "Encoder check " << (functor_multhpp::success ? "passed." : "failed.") << endl;
#if OUTPUT_DEBUGFILE
	debugfile << "Encoder check " << (functor_multhpp::success ? "passed." : "failed.") << endl;
#endif
#endif
	for (unsigned n = 0; n < N*Z; n++)
		// Perform BPSK and AWGN addition
		my[n] = (mx[n] ? -1 : 1) + grand()*sigma;	// 1->-1 and 0->1
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

	// Determine v(1)
	for (unsigned mi = 0; mi < Z; mi++)
	{
		functor_summsy::sum = Hp.pshift(0,0,mp,mi);	// P(i,k)v(0)
		Hs.iterY<functor_summsy>(mi);				// sigma P(i,j)u(j)
		mp[Z+mi] = functor_summsy::sum;				// p(1)
	}

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
}


unsigned LDPC::decode()
{
	// Set initial state
	for (unsigned n = 0; n < N*Z; n++)
	{
		functor_setq::l = 2.0/sigma/sigma*my[n]; //(col) output after AWGN and noise channel o/p
		//functor_setq::l = my[n];
		ml0[n] = functor_setq::l;	//ln^(0)
		ml[n] = functor_setq::l;	//ln (LLR update)
		H.iterX<functor_setq>(n);
	}

	// Iterative decoding
	for (unsigned i = 0; ; )
	{
#if OUTPUT_DEBUGFILE
		debugfile << "Before iteration " << i << ":" << endl;

		debugfile << "l:" << endl;
		outputLargeMatrix1<N*Z>(ml);

		debugfile << "q:" << endl;
		outputLargeMatrix2<N*Z,M*Z>(mq);

		debugfile.flush();

		if (i >= 1)
			return 0;
#endif
		// Update mr
		// rupdate_bp();
		rupdate_offms();

#if OUTPUT_DEBUGFILE
		debugfile << "After iteration " << i << ":" << endl;

		debugfile << "r:" << endl;
		outputLargeMatrix2<M*Z,N*Z>(mr);

		debugfile.flush();
#endif
		// Update mq and ml
		for (unsigned n = 0; n < Z*N; n++)
		{
			functor_sigmar::sigma = 0;

			H.iterX<functor_sigmar>(n);
			
			for (unsigned m = 0; m < Z*M; m++)
			{
				if (H.at(m, n))
				{
					mq[n][m] = mq0[n][m] + functor_sigmar::sigma; //variable node update 
					mq[n][m] -= mr[m][n]; //performs exlusion of m
				}
			}
			
			ml[n] = ml0[n] + functor_sigmar::sigma; //LLR update 

			mxhat[n] = ml[n] < 0; //hard decision
		}

		// Check that the decoding succeeded
		functor_multhxhat::nerrs = 0;
		H.multCol<functor_multhxhat>(mxhat);
		
		unsigned diff = 0;        //errors
		for (int j = 0; j < Z*K; j++)
			diff += mxhat[j] != ms[j];

#ifdef _DEBUG
		cout << "Iteration " << setw(3) << i << ": "
			<< setw(4) << functor_multhxhat::nerrs << " H*xhat errors, "
			<< setw(4) << diff << " x==xhat errors."
			<< endl;
#if OUTPUT_DEBUGFILE
		debugfile << "Iteration " << setw(3) << i << ": "
			<< setw(4) << functor_multhxhat::nerrs << " H*xhat errors, "
			<< setw(4) << diff << " x==xhat errors."
			<< endl;
#endif
#endif

		if (!functor_multhxhat::nerrs)
			return 0;

		if (++i > imax)
		{
			cout << "Maximum iteration reached; error." << endl;
#if OUTPUT_DEBUGFILE
			debugfile << "Maximum iteration reached; error." << endl;
#endif
			return 1;
		}
	}
}

void LDPC::rupdate_bp()
{
	// Update mr
	for (unsigned m = 0; m < Z*M; m++) //m outer dimension 
	{
		// Do the graph iteration to calculate the pi term without exclusion
		functor_r_bp::pi = 1;

		H.iterY<functor_r_bp>(m);

		for (unsigned n = 0; n < Z*N; n++) // n inner dimension
		{
			if (H.at(m, n))
			{
				long double pir = functor_r_bp::pi;
				const long double tanhr = tanh(mq[n][m]/2);
				if (tanhr)
					pir /= tanhr;
				else
				{
					n=n; // divide by 0 case
				}
				if (pir == 1)
				{
					mr[m][n] = numeric_limits<long double>::max(); // Explode!
				}
				else
				{
					const long double lnarg = (1+pir)/(1-pir);
					if (lnarg < 0)
					{
						n = n;	// Explode! -ve case
					}
					else
						mr[m][n] = log(lnarg);
				}
			}
			else
				mr[m][n] = 0;
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
		functor_r_offms::min0 = numeric_limits<long double>::max(); //use long double value as large value reference 
		functor_r_offms::min1 = functor_r_offms::min0;

		// Do the graph iteration to calculate the pi and min
		// terms without exclusion
		H.iterY<functor_r_offms>(m); //includes '1' ingnores
		
		for (int n = 0; n < Z*N; n++)
		{
			if (H.at(m, n))
			{
				long double pi = functor_r_offms::pi;

				// Perform exclusion
				long double qv = mq[n][m]; //Q matrix
				pi /= (qv > 0) ? 1 : -1; //performs exlusion
				long double qvmin = (fabs(qv) == functor_r_offms::min0) ?
					functor_r_offms::min1 : functor_r_offms::min0; //determines "true" the min

				mr[m][n] = pi * max((long double)0.0, qvmin - beta);//offset min sum calculation for r^(i)_(m,n)
			}
			else
				mr[m][n] = 0;
		}
	}
}