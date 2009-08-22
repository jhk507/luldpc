/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include <fstream>
#include <pthread.h>


#include "preachingbased.hpp"
#include "ldpcstate.hpp"

class LDPC
{
public:
	// The orthagonality error and message error histograms.
	OrthHistType orthhist[IMAX];
	MessHistType messhist[IMAX];

	// The performance histogram.
	PerfHistType perfhist;

	LDPCstate states[NTHREADS];

	int block;
	int nerrs;
	
	int totaliter[IMAX];	// Array to calculate average iteration
	
	pthread_mutex_t mutexcout;

#if OUTPUT_DEBUGFILE
	std::ofstream debugfile;
#endif

///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

	// Run the simulation
	void execute();

	void threadblock(LDPCstate *state);
	
	struct threadcontext
	{
		LDPC *ldpc;
		LDPCstate *state;
		pthread_t thread;
	};
	static void *threadproc(void *arg);
};

