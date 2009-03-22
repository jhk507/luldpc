/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include "automatrix.hpp"
#include "preachingbased.hpp"
#include "mtrand/MTRand_gaussian.hpp"

// Presets for half-rate
#define M 12	// Height of the unexpanded Preaching matrix
#define N 24	// Width of the unexpanded Preaching matrix

namespace LDPC
{
	// Initialize the simulation parameters
	void init();

	// Run the simulation
	void execute();

	// Computer the output of the encoder
	void encode();

	// Compute the output of the decoder
	// Returns true if no error
	bool decode();

	void rupdate_bp();

	void rupdate_offms();

	// Set the parity matrix based on the message
	void setParity();
}
