/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

// Presets for half-rate
#define M 12	// Height of the unexpanded Preaching matrix
#define N 24	// Width of the unexpanded Preaching matrix

namespace LDPC
{
	// Initialize the simulation parameters
	void setSnrDB(long double snrdbInit);

	// Run the simulation
	void execute();

	// Computer the output of the encoder
	void encode();

	// Set the parity matrix based on the message
	void setParity();

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
