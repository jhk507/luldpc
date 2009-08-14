/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include <fstream>

#include "preachingbased.hpp"
#include "ldpcstate.hpp"

class LDPC
{
public:
	// The orthagonality error and message error histograms.
	OrthHistType orthhist[IMAX];
	MessHistType messhist[IMAX];

	// The performance histogram.
	Histogram<NPERFBUCKETS, MAXPERFTIME> perfhist;

	LDPCstate states[NTHREADS];

	int block;
	int nerrs;

#if OUTPUT_DEBUGFILE
	std::ofstream debugfile;
#endif

///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

	// Run the simulation
	void execute();

	void threadblock(LDPCstate *state);
};
