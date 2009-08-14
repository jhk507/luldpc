/*
* $URL$
* $Date$
* $Rev$
*/

#include <iostream>

#include "itime.hpp"
#include "ldpc.hpp"
#include "histset.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

void LDPC::execute()
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

	for (DecodeMethod::Enum method = DecodeMethod::firstMethod; method < DecodeMethod::ndecodes; method++)
	{
		orthsets[method].init(decodeNames[method], "orth", orthhist);
		messsets[method].init(decodeNames[method], "mess", messhist);
	}

	ofstream perffile("hist_perf.tsv");

//	ITime runtimer;


	// The decode method loop
	for (DecodeMethod::Enum method = DecodeMethod::firstMethod;
		method < DecodeMethod::ndecodes; method++)
	{	
		if (method == DecodeMethod::bp || method == DecodeMethod::ms)
			continue; 
		cout << "Decoding using " << decodeNames[method] << " method...\n";

		// The SNR loop
		for (int snrindex = 0; snrindex < NSNRS; snrindex++)
		{
			block = 1;
			nerrs = 0;
			for (int ithread = 0; ithread < NTHREADS; ithread++)
			{
				LDPCstate *const state = states+ithread;
				state->init(method, snrindex);
				if (!ithread)
					cout << "\nSetting the SNR to " << state->snrdb << " dB...\n";
				threadblock(state);
			}
			cout << '\n';

			orthsets[method].writeLine(snrindex);
			messsets[method].writeLine(snrindex);

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
	for (DecodeMethod::Enum method = DecodeMethod::firstMethod; method < DecodeMethod::ndecodes; method++)
	{
		if (method == DecodeMethod::bp || method == DecodeMethod::ms)
			continue; 	
		decodeaxis << decodeNames[method] << '\n';
	}
	decodeaxis.close();

	ofstream snraxis("axis_snr.tsv");
	for (int snrindex = 0; snrindex < NSNRS; snrindex++)
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

void LDPC::threadblock(LDPCstate *state)
{
//	ITime blocktimer;
//	const double endtime = (RUNTIME-runtimer.get())/(NSNRS*(DecodeMethod::ndecodes-method)-snrindex);

	// The block loop
	for (; nerrs < NERRS; block++)
	{
		// Encode
		state->encode();

#if OUTPUT_DEBUGFILE
		debugfile << "Message:" << endl;
		outputLargeContiguous<K,Z>(ms, debugfile);

		debugfile << "Encoded parity bits:" << endl;
		outputLargeContiguous<M,Z>(mp, debugfile);
#endif

		// Decode
		if (!state->decode(orthhist, messhist, perfhist))
			nerrs++;

		if (!(block%10))
		{
			cout << "Block " << block
				<< "\tBLER=" << nerrs*100.0/block << "%        \r";
			cout.flush();
		}
	}
}
