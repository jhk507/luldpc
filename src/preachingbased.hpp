/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include "preaching.hpp"



//////////////////////////////////////////////////////////////////////////////

template <typename Elm, int Y, int X, int YRHO, int XRHO>
class PreachingBased
{
public:
	struct LinkedElm
	{
		Elm val;
		LinkedElm *nexty;
	//	inline operator Elm&() { return val; }
	};

public:
	PreachingBased(const Preaching<Y,X,YRHO,XRHO> &preachingInit);

	~PreachingBased();

	void output(std::ostream &out) const;

	template <typename Functor>
	inline void iterX(int x) const;

	template <typename Functor>
	inline void iterY(int y) const;

	template <typename Functor>
	inline void iterX2(int x, PreachingBased &p2) const;

	template <typename Functor>
	inline void iterY2(int y, PreachingBased &p2) const;

	template <typename Functor>
	inline void iterX3(int x, PreachingBased &p2, PreachingBased &p3) const;

	template <typename Functor>
	inline void iterY3(int y, PreachingBased &p2, PreachingBased &p3) const;

public:
	LinkedElm *data;
	LinkedElm *Vyc[Z*Y+1];
	LinkedElm *Vxc[Z*X];
	const Preaching<Y,X,YRHO,XRHO> &preaching;
};


template <typename Elm, int Y, int X, int YRHO, int XRHO>
PreachingBased<Elm,Y,X,YRHO,XRHO>::PreachingBased(const Preaching<Y,X,YRHO,XRHO> &preachingInit) :
	preaching(preachingInit)
{
	data = new LinkedElm[preaching.ones];

	static LinkedElm **rVyc, **rVxc;
	rVyc = Vyc;
	rVxc = Vxc;
	static const int (*rHyc)[Z*Y][YRHO+1];
	rHyc = &preaching.Hyc;

	LinkedElm *dati = data;
	for (int y = 0;; y++)
	{
		rVyc[y] = dati;
		if (y >= Z*Y)
			break;
		const int *pHyc = (*rHyc)[y];
		while (*pHyc++ >= 0)
			dati++;
	}

	for (int x = 0; x < Z*X; x++)
	{
		static LinkedElm *prev;
		prev = 0;
		struct functor_setlinks {
			static inline void callbackX(int y, int x) {
				LinkedElm *currentelm = rVyc[y];
				const int *currenth = (*rHyc)[y];
				while (*currenth != x)
				{
					currenth++;
					currentelm++;
				}
				if (prev)
					prev->nexty = currentelm;
				else
					rVxc[x] = currentelm;
				prev = currentelm;
			}
		};
		preaching.iterX<functor_setlinks>(x);
		prev->nexty = 0;
	}
}

template <typename Elm, int Y, int X, int YRHO, int XRHO>
PreachingBased<Elm,Y,X,YRHO,XRHO>::~PreachingBased()
{
	delete [] data;
}



template <typename Elm, int Y, int X, int YRHO, int XRHO>
template <typename Functor>
inline void PreachingBased<Elm,Y,X,YRHO,XRHO>::iterX(int x) const
{
	LinkedElm *e = Vxc[x];
	do
	{
		Functor::callback(e->val);
		e = e->nexty;
	} while (e);
}

template <typename Elm, int Y, int X, int YRHO, int XRHO>
template <typename Functor>
inline void PreachingBased<Elm,Y,X,YRHO,XRHO>::iterY(int y) const
{
	LinkedElm *e = Vyc[y];
	const LinkedElm *const end = Vyc[++y];
	do
	{
		Functor::callback(e->val);
		e++;
	} while (e != end);
}

template <typename Elm, int Y, int X, int YRHO, int XRHO>
template <typename Functor>
inline void PreachingBased<Elm,Y,X,YRHO,XRHO>::iterX2(int x, PreachingBased &p2) const
{
	LinkedElm *e1 = Vxc[x], *e2 = p2.Vxc[x];
	do
	{
		Functor::callback(e1->val, e2->val);
		e1 = e1->nexty;
		e2 = e2->nexty;
	} while (e1);
}

template <typename Elm, int Y, int X, int YRHO, int XRHO>
template <typename Functor>
inline void PreachingBased<Elm,Y,X,YRHO,XRHO>::iterY2(int y, PreachingBased &p2) const
{
	LinkedElm *e1 = Vyc[y], *e2 = p2.Vyc[y];
	const LinkedElm *const end = Vyc[++y];
	do
	{
		Functor::callback(e1->val, e2->val);
		e1++;
		e2++;
	} while (e1 != end);
}

template <typename Elm, int Y, int X, int YRHO, int XRHO>
template <typename Functor>
inline void PreachingBased<Elm,Y,X,YRHO,XRHO>::iterX3(int x, PreachingBased &p2, PreachingBased &p3) const
{
	LinkedElm *e1 = Vxc[x], *e2 = p2.Vxc[x], *e3 = p3.Vxc[x];
	do
	{
		Functor::callback(e1->val, e2->val, e3->val);
		e1 = e1->nexty;
		e2 = e2->nexty;
		e3 = e3->nexty;
	} while (e1);
}

template <typename Elm, int Y, int X, int YRHO, int XRHO>
template <typename Functor>
inline void PreachingBased<Elm,Y,X,YRHO,XRHO>::iterY3(int y, PreachingBased &p2, PreachingBased &p3) const
{
	LinkedElm *e1 = Vyc[y], *e2 = p2.Vyc[y];
	const LinkedElm *const end = Vyc[++y];
	do
	{
		Functor::callback(e1->val, e2->val, e3->val);
		e1++;
		e2++;
		e3++;
	} while (e1 != end);
}

template <typename Elm, int Y, int X, int YRHO, int XRHO>
void PreachingBased<Elm,Y,X,YRHO,XRHO>::output(std::ostream &out) const
{
	const LinkedElm *dati = data;
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
			out << std::endl;
	}
	out << std::endl;
}