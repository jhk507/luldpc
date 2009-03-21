#pragma once

#include "preaching.hpp"

//////////////////////////////////////////////////////////////////////////////

template <typename Elm, unsigned D>
class PreachingBased
{
protected:
	template <unsigned Y, unsigned X, unsigned RHO>
	PreachingBased(const Preaching<Y,X,RHO> &preaching)
	{
		data = new Elm[preaching.ones];
	}
public:
	~PreachingBased()
	{
		delete [] data;
	}

	Elm *Vc[Z*D];
	Elm *data;
};

//////////////////////////////////////////////////////////////////////////////

template <typename Elm, unsigned Y, unsigned X>
class PreachingBasedR :
	public PreachingBased<Elm, X>
{
public:
	template <unsigned RHO>
	PreachingBasedR(const Preaching<Y,X,RHO> &preaching) :
		PreachingBased(preaching)
	{
		Elm *dati = data;
		for (unsigned x = 0; x < Z*X; x++)
		{
			Vc[x] = dati;
			int *pHxc = preaching.Hxc[x];
			while (*pHxc++ >= 0)
				dati++;
		}
	}
};

//////////////////////////////////////////////////////////////////////////////


template <typename Elm, unsigned Y, unsigned X>
class PreachingBasedQ :
	public PreachingBased<Elm, Y>
{
public:
	template <unsigned RHO>
	PreachingBasedQ(const Preaching<Y,X,RHO> &preaching) :
		PreachingBased(preaching)
	{
		Elm *dati = data;
		for (unsigned y = 0; y < Z*Y; y++)
		{
			Vc[y] = dati;
			int *pHyc = preaching.Hyc[x];
			while (*pHyc++ >= 0)
				dati++;
		}
	}
};
