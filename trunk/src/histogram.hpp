/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include <ostream>

template <int nBuckets, int valMax, int valSection, int nTrials>
class Histogram
{
public:
	inline void report(int val)
	{
		int bucket = val ? 1+(nBuckets-1)*(val-1)/valSection: 0;
		if (bucket >= nBuckets)
			bucket = nBuckets-1;
		buckets[bucket]++;
	}

	static void outputHeader(std::ostream &out)
	{
		out << '0';
		for (int b = 1; b < nBuckets; b++)
		{
			const long double v =
				(valSection*(b-1)/(long double)(nBuckets-1) + 1)/valMax;
			out << '\t' << v;
		}
		out << '\n';
	}

	void output(std::ostream &out) const
	{
		for (int b = 0; b < nBuckets; b++)
			out << buckets[b]/(long double)nTrials << '\t';
		out << '\n';
	}

private:
	int buckets[nBuckets];
};
