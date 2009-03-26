/*
* $URL$
* $Date$
* $Rev$
*/

#include <iostream>
#if defined(_MSC_VER)
#include <windows.h>
#endif

#include "ldpc.hpp"

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

using namespace std;

int main()
{
#if defined(_MSC_VER) && defined(_DEBUG)
	// Enable floating-point exceptions for debugging purposes
	enableFPEs();
#endif

#if defined(WIN32) && !defined(_DEBUG)
	SetPriorityClass(GetCurrentProcess(), ABOVE_NORMAL_PRIORITY_CLASS);
	SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_ABOVE_NORMAL);

	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);
	LARGE_INTEGER start;
	QueryPerformanceCounter(&start);
#endif

	LDPC::execute();

#if defined(WIN32) && !defined(_DEBUG)
	LARGE_INTEGER end;
	QueryPerformanceCounter(&end);
	end.QuadPart -= start.QuadPart;
	cout << "\n\nRuntime: " << ((long double)end.QuadPart)/freq.QuadPart << "s\n";
#endif

	return 0;
}
