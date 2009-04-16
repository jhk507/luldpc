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

#define RUNTIME      60		// The number of seconds to run (approximate)
#define NERRBUCKETS  25		// The number of error histogram buckets
#define NPERFBUCKETS 100	// The number of performance histogram buckets
#define MAXPERFTIME  2000	// The maximum performance histogram duration, in us

// The orthagonality error and message error histograms.
// The template parameters are the number of histogram buckets, the full size
// of the data range, and the desired portion of the data range to examine.
typedef ValNormalizedHistogram<NERRBUCKETS, M*Z, (int)(M*Z*0.33)> OrthHistType;
typedef ValNormalizedHistogram<NERRBUCKETS, N*Z, (int)(N*Z*0.06)> MessHistType;
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
