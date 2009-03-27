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
		buckets[getBucket(val)]++;
	}

	static inline double getNormalizedValFloor(int b)
	{
		if (b)
			return (valSection*(b-1)/(double)(nBuckets-1) + 1)/valMax;
		return 0;
	}

	static inline int getBucket(int val)
	{
		if (!val)
			return 0;
		const int bucket = 1+(nBuckets-1)*(val-1)/valSection;
		if (bucket < nBuckets)
			return bucket;
		return nBuckets-1;
	}

	inline double getNormalizedFreq(int b) const
	{
		return buckets[b]/(double)nTrials;
	}

private:
	int buckets[nBuckets];
};
