/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

namespace LDPC
{
///////////////////////////////////////////////////////////////////////////////
// Globals ////////////////////////////////////////////////////////////////////

// Presets for half-rate
#define M 12		// Height of the unexpanded Preaching matrix
#define N 24		// Width of the unexpanded Preaching matrix

#define IMAX 50		// The maximum number of decode iterations

// Matrix sparsity parameters
#define RHO_H_Y  7
#define RHO_H_X  6
#define RHO_HS_Y 5
#define RHO_HS_X 6
#define RHO_HP_Y 3
#define RHO_HP_X 3

// Data for unexpanded half-rate Preaching matrix H
const extern int Ha[M][N];

const extern Preaching<M,N,RHO_H_Y, RHO_H_X>  H;	// Unexpanded half-rate Preaching matrix H
const extern Preaching<M,K,RHO_HS_Y,RHO_HS_X> Hs;	// Half-rate Preaching matrix H (first half)
const extern Preaching<M,M,RHO_HP_Y,RHO_HP_X> Hp;	// Half-rate Preaching matrix H (second half, for parity)

extern PreachingBased<double, M,N,RHO_H_Y,RHO_H_X> mr;	// R matrix
extern PreachingBased<double, M,N,RHO_H_Y,RHO_H_X> mq;	// Q matrix

// The decode method.
namespace DecodeMethod
{
enum Enum
{
	firstMethod = 0,	// (The first method available)
	ms = 0,				// Min sum
	ms_sc,
	offms,				// Offset min sum
	offms_sc,
	nms,				// Normalized min sum
	nms_sc,
	bp,					// Belief propagation
	ndecodes			// (The number of decoding algorithms)
};
}
extern DecodeMethod::Enum method;

extern const char *const decodeNames[DecodeMethod::ndecodes];

///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

inline void operator++(DecodeMethod::Enum &incmethod, int)
{
	incmethod = (DecodeMethod::Enum)((int)incmethod + 1);
}

// Compute the output of the decoder
// Returns true if no error
bool decode();

// Set the initial decoder state
void decode_initial();

// Update the R matrix (belief propagation method)
void rupdate_bp();

// Update the R matrix (minsum method)
template <DecodeMethod::Enum msMethod, bool sc>
void rupdate_ms();

// Update the Q and L matrices
template <bool sc>
void qlupdate();

}
