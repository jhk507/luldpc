#include <iostream>

#include "ldpc.hpp"

using namespace std;

#if defined(_MSC_VER) && defined(_DEBUG)

#include <float.h>

void enableFPEs()
{
	_clearfp(); //Always call _clearfp before setting the control
				//word
	//Because the second parameter in the following call is 0, it
	//only returns the floating-point control word
	unsigned int cw = _controlfp(0, 0); //Get the default control
										//word
	//Set the exception masks off for exceptions that you want to
	//trap.  When a mask bit is set, the corresponding floating-point
	//exception is //blocked from being generating.
	cw &=~(EM_OVERFLOW // |EM_UNDERFLOW
		| EM_ZERODIVIDE|EM_DENORMAL|EM_INVALID);
	//For any bit in the second parameter (mask) that is 1, the
	//corresponding bit in the first parameter is used to update
	//the control word.
	unsigned int cwOriginal = _controlfp(cw, MCW_EM); //Set it.
								//MCW_EM is defined in float.h.
								//Restore the original value when done:
								//_controlfp(cwOriginal, MCW_EM);

	unsigned long cntrReg;
	__asm stmxcsr [cntrReg]       //Get MXCSR register
	cntrReg &= 0xFFFFFF7F & ~(2<<11); //bit 7 - invalid instruction mask
								  //bit 9  - divide-by-zero mask
								  //bit 10 - overflow mask
								  //bit 11 - underflow mask
	__asm ldmxcsr [cntrReg]       //Load MXCSR register
}
#endif


int main()
{
#if defined(_MSC_VER) && defined(_DEBUG)
	// Enable floating-point exceptions for debugging purposes
	enableFPEs();
#endif

	LDPC::execute();

	return 0;
}
