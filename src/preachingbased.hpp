#pragma once

#include <ostream>

#include "preaching.hpp"

//////////////////////////////////////////////////////////////////////////////

template <typename Elm, int Y, int X, int RHO>
class PreachingBased
{
public:
	PreachingBased(const Preaching<Y,X,RHO> &preachingInit) :
		preaching(preachingInit)
	{
		data = new Elm[preaching.ones];

		Elm *dati = data;
		for (int y = 0; y < Z*Y; y++)
		{
			Vyc[y] = dati;
			const int *pHyc = preaching.Hyc[y];
			while (*pHyc++ >= 0)
				dati++;
		}

		dati = data;
		for (int x = 0; x < Z*X; x++)
		{
			Vxc[x] = dati;
			const int *pHxc = preaching.Hxc[x];
			while (*pHxc++ >= 0)
				dati++;
		}
	}

	~PreachingBased()
	{
		delete [] data;
	}

	void output(std::ostream &out) const
	{
		const Elm *dati = data;
		for (int y = 0; y < Z*Y; y++)
		{
			const int *pHyc = preaching.Hyc[y];
			int xh = *pHyc;
			for (int x = 0; x < Z*X; x++)
			{
				if (xh == x)
				{
					xh = *++pHyc;
					out << *++dati << '\t';
				}
				else
					out << "0\t";

				if (x%Z == Z-1)
					out << '\t';
			}
			out << std::endl;
			if (y%Z == Z-1)
				out << endl;
		}
		out << std::endl;
	}

	Elm *data;
	Elm *Vyc[Z*Y];
	Elm *Vxc[Z*X];
	const Preaching<Y,X,RHO> &preaching;
};
