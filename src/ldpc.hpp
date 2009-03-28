/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include "histogram.hpp"
#include "preachingbased.hpp"
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

const extern Preaching<M,N,RHO_H_Y, RHO_H_X>  H;	// Unexpanded half-rate Preaching matrix H
const extern Preaching<M,K,RHO_HS_Y,RHO_HS_X> Hs;	// Half-rate Preaching matrix H (first half)
const extern Preaching<M,M,RHO_HP_Y,RHO_HP_X> Hp;	// Half-rate Preaching matrix H (second half, for parity)

extern PreachingBased<double, M,N,RHO_H_Y,RHO_H_X> mr;	// R matrix
extern PreachingBased<double, M,N,RHO_H_Y,RHO_H_X> mq;	// Q matrix

// The orthagonality error and message error histograms.
// The template parameters are the number of histogram buckets, the full size
// of the data range, and the desired portion of the data range to examine.
typedef Histogram<NBUCKETS, M*Z, (int)(M*Z*0.33)> OrthHistType;
typedef Histogram<NBUCKETS, N*Z, (int)(N*Z*0.06)> MessHistType;
extern OrthHistType orthhist[IMAX];
extern MessHistType messhist[IMAX];

///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

// Initialize the simulation parameters
void setSnrDB(double snrdbInit);

// Run the simulation
void execute();

}
