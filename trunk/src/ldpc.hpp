/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include "histogram.hpp"
#include "preaching.hpp"
#include "decode.hpp"

namespace LDPC
{
///////////////////////////////////////////////////////////////////////////////
// Globals ////////////////////////////////////////////////////////////////////

// Presets for half-rate
#define M 12	// Height of the unexpanded Preaching matrix
#define N 24	// Width of the unexpanded Preaching matrix

// Matrix sparsity parameters
#define RHO_H_Y  7
#define RHO_H_X  6
#define RHO_HS_Y 5
#define RHO_HS_X 6
#define RHO_HP_Y 3
#define RHO_HP_X 3

#define OUTPUT_DEBUGFILE 0	// Enable to output data to a debug file

#define NBUCKETS 25	// The number of histogram buckets

const extern Preaching<M,N,RHO_H_Y, RHO_H_X> H;	// Unexpanded half-rate Preaching matrix H

extern bool mx[N*Z];	// (col) Combination of ms and mp

// Set the aliases into mx
extern bool (&ms)[K*Z];	// (col) Message
extern bool (&mp)[M*Z];	// (col) Generated parity

// The orthagonality error and message error histograms.
// The template parameters are the number of histogram buckets, the full size
// of the data range, and the desired portion of the data range to examine.
extern Histogram<NBUCKETS, M*Z, (int)(M*Z*0.33)> orthhist[IMAX];
extern Histogram<NBUCKETS, N*Z, (int)(N*Z*0.06)> messhist[IMAX];

///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

// Initialize the simulation parameters
void setSnrDB(double snrdbInit);

// Run the simulation
void execute();

}
