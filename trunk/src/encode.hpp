/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include "ldpc.hpp"

namespace LDPC
{
///////////////////////////////////////////////////////////////////////////////
// Globals ////////////////////////////////////////////////////////////////////

#define NSNRS (sizeof(snrs)/sizeof(*snrs))	// The number of SNRs to try
#define DEFAULTSNR 5						// The index of the default SNR.

extern double sigma;

extern double my[N*Z];		// (col) Encoder output after AWGN

// The SNRs to try.
extern const double snrs[] = { 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6 };

extern int snrindex;								// The current SNR index


///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

// Computer the output of the encoder
void encode();

// Set the parity matrix based on the message
void setParity();

}
