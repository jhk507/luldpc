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
	ITime()
	{
		reset();
	}

	void reset();

	inline double get() const
	{
#ifdef WIN32
	LARGE_INTEGER end;
	QueryPerformanceCounter(&end);
	end.QuadPart -= start.QuadPart;
	return end.QuadPart/(double)freq.QuadPart;
#else
	timeval end;
	gettimeofday(&end, 0);
	return end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)/1000000.0;
#endif
	}

private:
#ifdef WIN32
	static bool freqInitialized;
	static LARGE_INTEGER freq;

	LARGE_INTEGER start;
#else
	timeval start;
#endif
};

class Profiler : public ITime
{
public:
	void start()
	{
		reset();
	}

	void stop()
	{
		const double g = get();
		const double t = average*count;
		count++;
		average = (g+t)/count;
	};

private:
	double average;
	int count;
};

}
