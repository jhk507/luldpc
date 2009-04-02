/*
* $URL$
* $Date$
* $Rev$
*/

#ifdef WIN32
#include <windows.h>
#undef max
#else
#include <sys/time.h>
#endif

#pragma once

namespace LDPC
{

class ITime
{
public:
	ITime();

	double get() const;

private:
#ifdef WIN32
	static bool freqInitialized;
	static LARGE_INTEGER freq;

	LARGE_INTEGER start;
#else
	timeval start;
#endif
};

}
