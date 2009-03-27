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

#define IMAX 50		// The maximum number of decode iterations

///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

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
