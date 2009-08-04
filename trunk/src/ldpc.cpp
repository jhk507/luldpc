/*
* $URL$
* $Date$
* $Rev$

Warning: This code is VERY NON-REENTRANT. This is deliberate to facilitate
access to variables from functors and to greatly increase performance by
cutting down on frame pointer generation. In short, no multithreading allowed.
*/

#include <iostream>

#include "itime.hpp"
#include "encode.hpp"
#include "histset.hpp"

using namespace std;

namespace LDPC
{

///////////////////////////////////////////////////////////////////////////////
// Globals ////////////////////////////////////////////////////////////////////

// The orthagonality error and message error histograms.
OrthHistType orthhist[IMAX];
MessHistType messhist[IMAX];

// The performance histogram.
Histogram<NPERFBUCKETS, MAXPERFTIME> perfhist;


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

	HistogramSet<OrthHistType> orthsets[DecodeMethod::ndecodes];
	HistogramSet<MessHistType> messsets[DecodeMethod::ndecodes];

	for (method = DecodeMethod::firstMethod; method < DecodeMethod::ndecodes; method++)
	{
		orthsets[method].init(decodeNames[method], "orth", orthhist);
		messsets[method].init(decodeNames[method], "mess", messhist);
	}

	ofstream perffile("hist_perf.tsv");

//	ITime runtimer;

	// The decode method loop
	for (method = DecodeMethod::firstMethod; method < DecodeMethod::ndecodes; method++)
	//method = DecodeMethod::ms;
	{
		if (method == DecodeMethod::bp)
			continue;
		if (method == DecodeMethod::ms)
			continue; 

		cout << "Decoding using " << decodeNames[method] << " method...\n";

		// The SNR loop
		for (snrindex = 0; snrindex < NSNRS; snrindex++)
		{
			// Set the signal-to-noise ratio
			setSnrDB(snrs[snrindex]);

			int nerrs = 0;	// The number of block errors

//			ITime blocktimer;
//			const double endtime = (RUNTIME-runtimer.get())/(NSNRS*(DecodeMethod::ndecodes-method)-snrindex);

			// The block loop
			for (int b = 1; nerrs<NERRS ; b++)
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
					cout << "Block " << b
						<< "\tBLER=" << nerrs*100.0/b << "%        \r";
					cout.flush();
				}
							
			}
			cout<<'\n';

			orthsets[method].writeLine();
			messsets[method].writeLine();

			for (int i = 0; i < IMAX; i++)
			{
				orthhist[i].reset();
				messhist[i].reset();
			}
		}
		cout << '\n';

		for (int b = 0; b < NPERFBUCKETS; b++)
			perffile << perfhist.getNormalizedFreq(b) << '\t';
		perffile << '\n';
		perfhist.reset();
	}

	// Output the rest of the TSV data.
	cout << "Generating remaining data files...\n";

	ofstream decodeaxis("axis_decode.tsv");
	for (method = DecodeMethod::firstMethod; method < DecodeMethod::ndecodes; method++)
	{	if (method == DecodeMethod::bp)
			continue;
		if (method == DecodeMethod::ms)
			continue;   
	//	if (method == DecodeMethod::nms)
	//		continue; 
	//	if (method == DecodeMethod::nms_sc)
	//		continue; 
	//	if (method == DecodeMethod::offms)
	//		continue; 
	//	if (method == DecodeMethod::offms_sc)
	//		continue; 
	//	if (method == DecodeMethod::v_off_ms)
	//		continue; 
	
		decodeaxis << decodeNames[method] << '\n';
	}
		decodeaxis.close();

	ofstream snraxis("axis_snr.tsv");
	for (snrindex = 0; snrindex < NSNRS; snrindex++)
		snraxis << snrs[snrindex] << '\n';
	snraxis.close();

	ofstream iteraxis("axis_iter.tsv");
	for (int i = 1; i <= IMAX; i++)
		iteraxis << i << '\n';
	iteraxis.close();

	ofstream orthaxis("axis_err_orth.tsv");
	ofstream messaxis("axis_err_mess.tsv");
	for (int b = 0; b < NERRBUCKETS; b++)
	{
		orthaxis << orthhist->getValFloor(b) << '\n';
		messaxis << messhist->getValFloor(b) << '\n';
	}
	orthaxis.close();
	messaxis.close();

	ofstream perfaxis("axis_perf.tsv");
	for (int b = 0; b < NPERFBUCKETS; b++)
		perfaxis << perfhist.getValFloor(b) << '\n';
	perfaxis.close();
}

}
