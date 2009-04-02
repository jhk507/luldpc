/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

template <int nBuckets, int valMax>
class Histogram
{
public:
	Histogram()
	{
		reset();
	}

	void reset()
	{
		ntrials = 0;
		for (int b = 0; b < nBuckets; b++)
			buckets[b] = 0;
	}

	inline void report(int val)
	{
		buckets[getBucket(val)]++;
		ntrials++;
	}

	inline double getNormalizedFreq(int b) const
	{
		return buckets[b]/(double)ntrials;
	}

	virtual double getValFloor(int b) const
	{
		return valMax/0.999*b/nBuckets;
	}

	virtual int getBucket(int val) const
	{
		const int b = (int)(val*0.999*nBuckets/valMax);
		if (b < nBuckets)
			return b;
		return nBuckets-1;
	}

private:
	int buckets[nBuckets];
	int ntrials;
};


template <int nBuckets, int valMax, int valSection>
class ValNormalizedHistogram : public Histogram<nBuckets, valMax>
{
public:
	virtual double getValFloor(int b) const
	{
		if (b)
			return (valSection*(b-1)/(double)(nBuckets-1) + 1)/valMax;
		return 0;
	}

	virtual int getBucket(int val) const
	{
		if (!val)
			return 0;
		const int bucket = 1+(nBuckets-1)*(val-1)/valSection;
		if (bucket < nBuckets)
			return bucket;
		return nBuckets-1;
	}
};

