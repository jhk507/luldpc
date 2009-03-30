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

#define NSNRS 4			// The number of SNRs to try
#define DEFAULTSNR 2	// The index of the default SNR.

extern double sigma;

extern bool mx[N*Z];	// (col) Combination of ms and mp
extern double my[N*Z];	// (col) Encoder output after AWGN

// Set the aliases into mx
extern bool (&ms)[K*Z];	// (col) Message
extern bool (&mp)[M*Z];	// (col) Generated parity

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
