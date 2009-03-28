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

#define NSNRS 10		// The number of SNRs to try
#define DEFAULTSNR 5	// The index of the default SNR.

extern double sigma;

extern double my[N*Z];	// (col) Encoder output after AWGN

// The SNRs to try.
extern const double snrs[NSNRS];

extern int snrindex;	// The current SNR index


///////////////////////////////////////////////////////////////////////////////
// Functions //////////////////////////////////////////////////////////////////

// Computer the output of the encoder
void encode();

// Set the parity matrix based on the message
void setParity();

}
