/*
* $URL$
* $Date$
* $Rev$
*/

#include "itime.hpp"

namespace LDPC
{

#ifdef WIN32
LARGE_INTEGER ITime::freq;
bool ITime::freqInitialized = QueryPerformanceFrequency(&ITime::freq);

ITime::ITime()
{
	QueryPerformanceCounter(&start);
}

double ITime::get() const
{
	LARGE_INTEGER end;
	QueryPerformanceCounter(&end);
	end.QuadPart -= start.QuadPart;
	return end.QuadPart/(double)freq.QuadPart;
}

#else

ITime::ITime()
{
	gettimeofday(&start, 0);
}

double ITime::get() const
{
	timeval end;
	gettimeofday(&end, 0);
	return end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)/1000000.0;
}

#endif

};
