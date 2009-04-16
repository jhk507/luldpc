/*
* $URL$
* $Date$
* $Rev$
*/

#include <iostream>
#include <cmath>

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

	inline void reset()
	{
#ifdef WIN32
		QueryPerformanceCounter(&start);
#else
		gettimeofday(&start, 0);
#endif
	}

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

	void outputLong(std::ostream &out) const
	{
		const double now = get();
		const int hours = now/3600;
		const int mins = (int)(now/60) % 60;
		const double secs = fmod(now, 60);

		out << hours << ':';
		if (mins < 10)
			out << '0';
		out << mins << ':';
		if (secs < 10)
			out << '0';
		out << secs << '\n';
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
	inline void start()
	{
		reset();
	}

	inline void stop()
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
