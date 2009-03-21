#pragma once

#include "preaching.hpp"

//////////////////////////////////////////////////////////////////////////////

template <typename Elm>
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

	Elm *data;
};

//////////////////////////////////////////////////////////////////////////////

template <typename Elm, unsigned Y>
class PreachingBasedY :
	public PreachingBased<Elm>
{
public:
	template <unsigned X, unsigned RHO>
	PreachingBasedY(const Preaching<Y,X,RHO> &preaching) :
		PreachingBased<Elm>(preaching)
	{
		Elm *dati = PreachingBased<Elm>::data;
		for (unsigned y = 0; y < Z*Y; y++)
		{
			Vcy[y] = dati;
			int *pHyc = preaching.Hyc[y];
			while (*pHyc++ >= 0)
				dati++;
		}
	}

	Elm *Vcy[Z*Y];
};

//////////////////////////////////////////////////////////////////////////////

template <typename Elm, unsigned X>
class PreachingBasedX :
	virtual public PreachingBased<Elm>
{
public:
	template <unsigned Y, unsigned RHO>
	PreachingBasedX(const Preaching<Y,X,RHO> &preaching) :
		PreachingBased<Elm>(preaching)
	{
		Elm *dati = PreachingBased<Elm>::data;
		for (unsigned x = 0; x < Z*X; x++)
		{
			Vcx[x] = dati;
			int *pHxc = preaching.Hxc[x];
			while (*pHxc++ >= 0)
				dati++;
		}
	}

	Elm *Vcx[Z*X];
};

//////////////////////////////////////////////////////////////////////////////

template <typename Elm, unsigned Y, unsigned X>
class PreachingBasedR :
	public PreachingBasedX<Elm, X>
{
public:
	template <unsigned RHO>
	PreachingBasedR(const Preaching<Y,X,RHO> &preaching) :
		PreachingBasedX<Elm, X>(preaching)
	{
	}
};

template <typename Elm, unsigned Y, unsigned X>
class PreachingBasedQ :
	public PreachingBasedY<Elm, Y>,
	public PreachingBasedX<Elm, X>
{
public:
	template <unsigned RHO>
	PreachingBasedQ(const Preaching<Y,X,RHO> &preaching) :
		PreachingBasedY<Elm, Y>(preaching),
		PreachingBasedX<Elm, X>(preaching)
	{
	}
};
