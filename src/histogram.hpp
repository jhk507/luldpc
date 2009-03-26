/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

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

	static inline double getNormalizedBucket(int b)
	{
		if (b)
			return (valSection*(b-1)/(double)(nBuckets-1) + 1)/valMax;
		return 0;
	}

	inline double getNormalizedFreq(int b) const
	{
		return buckets[b]/(double)nTrials;
	}

private:
	int buckets[nBuckets];
};
