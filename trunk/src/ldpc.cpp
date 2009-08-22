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

	HistogramSet<OrthHistType> orthsets[DecodeMethod::ndecodes];
	HistogramSet<MessHistType> messsets[DecodeMethod::ndecodes];

	for (DecodeMethod::Enum method = DecodeMethod::firstMethod; method < DecodeMethod::ndecodes; method++)
	{
		orthsets[method].init(decodeNames[method], "orth", orthhist);
		messsets[method].init(decodeNames[method], "mess", messhist);
	}

	ofstream perffile("hist_perf.tsv");

//	ITime runtimer;

	threadcontext contexts[NTHREADS];
	for (int t = 0; t < NTHREADS; t++)
	{
		contexts[t].ldpc = this;
		contexts[t].state = &states[t];
	}
	
	if (pthread_mutex_init(&mutexcout, 0))
	{
		cerr << "Could not create console output mutex!\n";
		return;
	}

	// The decode method loop
	for (DecodeMethod::Enum method = DecodeMethod::firstMethod;
		method < DecodeMethod::ndecodes; method++)
	{
		if (method == DecodeMethod::bp || method == DecodeMethod::ms)
			continue; 
		cout << "Decoding using " << decodeNames[method] << " method...\n";
		
		string iterfilename = "hist_iter_";
		iterfilename += decodeNames[method];
		iterfilename += ".tsv";
		ofstream iterfile(iterfilename.c_str());
		iterfile << setprecision(10);

		// The SNR loop
		for (int snrindex = 0; snrindex < NSNRS; snrindex++)
		{
			block = 1;
			nerrs = 0;
			
			for (int i = 0; i < IMAX; i++)
				totaliter[i] = 0;
				
			for (int ithread = 0; ithread < NTHREADS; ithread++)
			{
				LDPCstate *const state = states+ithread;
				state->init(method, snrindex);
				if (!ithread)
					cout << "\nSetting the SNR to " << state->snrdb << " dB...\n";
				//threadblock(state);
				
				if (pthread_create(&contexts[ithread].thread, 0, threadproc, contexts+ithread))
				{
					cerr << "Thread creation failed!\n";
					return;
				}
			}
			
			for (int ithread = 0; ithread < NTHREADS; ithread++)
			{
				void *status;
				if (pthread_join(contexts[ithread].thread, &status))
				{
					cerr << "Could not join thread!\n";
					return;
				}
			}
			
			cout << '\n';

			orthsets[method].writeLine(snrindex);
			messsets[method].writeLine(snrindex);

			for (int i = 0; i < IMAX; i++)
			{
				orthhist[i].reset();
				messhist[i].reset();
				
				const long double aveiter = totaliter[i]/((long double)block);
				iterfile << aveiter << '\t';
			}
			iterfile << '\n';
		}
		cout << '\n';

		for (int b = 0; b < NPERFBUCKETS; b++)
			perffile << perfhist.getNormalizedFreq(b) << '\t';
		perffile << '\n';
		perfhist.reset();
	}

	if (pthread_mutex_destroy(&mutexcout))
		cerr << "Could not destroy console output mutex!\n";
	
	ofstream decodeaxis("axis_decode.tsv");
	for (DecodeMethod::Enum method = DecodeMethod::firstMethod; method < DecodeMethod::ndecodes; method++)
	{
		if (method == DecodeMethod::bp || method == DecodeMethod::ms)
			continue; 	
		decodeaxis << decodeNames[method] << '\n';
	}
	decodeaxis.close();

	ofstream perfaxis("axis_perf.tsv");
	for (int b = 0; b < NPERFBUCKETS; b++)
		perfaxis << perfhist.getValFloor(b) << '\n';
	perfaxis.close();
}

void *LDPC::threadproc(void *arg)
{
	threadcontext *const instance = (threadcontext*)arg;
	instance->ldpc->threadblock(instance->state);
	return 0;
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
			
		// Add to totaliter array for "average iteration" calculation
		int i;
		for (i = 0; i < state->iter; i++)
			totaliter[i] += i;
		for (; i < IMAX; i++)
			totaliter[i] += state->iter;

		if (!(block%100))
		{
			if (!pthread_mutex_trylock(&mutexcout))
			{
				cout << "Block " << block
					<< "    BLER=" << nerrs*100.0/block << "%    " << nerrs << " errors       \r";
				cout.flush();
				pthread_mutex_unlock(&mutexcout);
			}
		}
		
//		usleep(1000);
	}
}
