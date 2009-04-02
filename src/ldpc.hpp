/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include <fstream>

#include "histogram.hpp"
#include "preachingbased.hpp"
#include "decode.hpp"

namespace LDPC
{
///////////////////////////////////////////////////////////////////////////////
// Globals ////////////////////////////////////////////////////////////////////

#define OUTPUT_DEBUGFILE 0	// Enable to output data to a debug file

#define NBLOCKS 100			// The number of blocks to run
#define NERRBUCKETS 25		// The number of error histogram buckets
#define NPERFBUCKETS 100
#define MAXPERFTIME  1000

// The orthagonality error and message error histograms.
// The template parameters are the number of histogram buckets, the full size
// of the data range, and the desired portion of the data range to examine.
typedef ValNormalizedHistogram<NERRBUCKETS, M*Z, (M*Z*0.33)> OrthHistType;
typedef ValNormalizedHistogram<NERRBUCKETS, N*Z, (N*Z*0.06)> MessHistType;
extern OrthHistType orthhist[IMAX];
extern MessHistType messhist[IMAX];

// The performance histogram.
extern Histogram<NPERFBUCKETS, MAXPERFTIME> perfhist;

#if OUTPUT_DEBUGFILE
extern std::ofstream debugfile;
#endif

///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

// Initialize the simulation parameters
void setSnrDB(double snrdbInit);

// Run the simulation
void execute();

}
