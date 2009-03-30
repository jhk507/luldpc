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
#define M 12	// Height of the unexpanded Preaching matrix
#define N 24	// Width of the unexpanded Preaching matrix

#define IMAX 100	// The maximum number of decode iterations

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
enum DecodeMethod
{
	firstMethod = 0,	// (The first method available)
	bp = 0,				// Belief propagation
	offms,				// Offset min sum
	ndecodes			// (The number of decoding algorithms)
};
extern DecodeMethod method;

extern const char *const decodeNames[ndecodes];

///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

inline void operator++(DecodeMethod &incmethod, int)
{
	incmethod = (DecodeMethod)((int)incmethod + 1);
}

// Compute the output of the decoder
// Returns true if no error
bool decode();

// Set the initial decoder state
void decode_initial();

// Update the R matrix (belief propagation method)
void rupdate_bp();

// Update the R matrix (offset minsum method)
void rupdate_offms();

// Update the Q and L matrices
void qlupdate();

}
