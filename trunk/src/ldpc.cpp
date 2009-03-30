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

	for (method = firstMethod; method < ndecodes; method++)
	{
		orthsets[method].init(decodeNames[method], "orth", orthhist);
		messsets[method].init(decodeNames[method], "mess", messhist);
	}

	// The decode method loop
	for (method = firstMethod; method < ndecodes; method++)
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
						<< "\tBLER=" << nerrs*100.0/b << "%        \r";
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

	ofstream decodeaxis("axis_decode.tsv");
	for (method = firstMethod; method < ndecodes; method++)
		decodeaxis << decodeNames[method] << '\n';
	decodeaxis.close();

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
