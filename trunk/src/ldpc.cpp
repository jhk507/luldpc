/*
* $URL$
* $Date$
* $Rev$

Warning: This code is VERY NON-REENTRANT. This is deliberate to facilitate
access to variables from functors and to greatly increase performance by
cutting down on frame pointer generation. In short, no multithreading allowed.
*/

#include <iostream>

#include "encode.hpp"
#include "histset.hpp"

using namespace std;

namespace LDPC
{

///////////////////////////////////////////////////////////////////////////////
// Globals ////////////////////////////////////////////////////////////////////

#define NBLOCKS 100	// The number of blocks to run

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

// The orthagonality error and message error histograms.
// The template parameters are the number of histogram buckets, the full size
// of the data range, and the desired portion of the data range to examine.
OrthHistType orthhist[IMAX];
MessHistType messhist[IMAX];

#if OUTPUT_DEBUGFILE
ofstream debugfile("debugfile.tsv");
#endif

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

	HistogramSet<OrthHistType> orthsets[ndecodes];
	HistogramSet<MessHistType> messsets[ndecodes];

	int &imethod = (int&)method;
	for (method = firstMethod; method < ndecodes; imethod++)
	{
		orthsets[method].init(decodeNames[method], "orth", orthhist);
		messsets[method].init(decodeNames[method], "mess", messhist);
	}

	// The decode method loop
	for (method = firstMethod; method < ndecodes; imethod++)
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
				{
					cout << "Block " << b << '/' << NBLOCKS
						<< "\tBLER=" << nerrs*100.0/NBLOCKS << "%    \r";
					cout.flush();
				}
			}
			cout << '\n';

			orthsets[method].writeLine();
			messsets[method].writeLine();

			for (int i = 0; i < IMAX; i++)
			{
				orthhist[i].reset();
				messhist[i].reset();
			}
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
		orthaxis << orthhist->getNormalizedValFloor(b) << '\n';
		messaxis << messhist->getNormalizedValFloor(b) << '\n';
	}
	orthaxis.close();
	messaxis.close();
}

}
