/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include "preaching.hpp"



//////////////////////////////////////////////////////////////////////////////

// A matrix with the same dimensions as the Preaching matrix, and with values
// in the same places of the stored sparse matrix as in those of the Preaching
// matrix.
template <typename Elm, int Y, int X, int YRHO, int XRHO>
class PreachingBased
{
public:
	// An element structure containing the value at certain coordinates
	// as well as a pointer to the next value in the y direction.
	struct LinkedElm
	{
		Elm val;
		LinkedElm *nexty;
	};

public:
	// The constructor. Takes a reference to an existing Preaching object.
	PreachingBased(const Preaching<Y,X,YRHO,XRHO> &preachingInit);

	// The destructor.
	~PreachingBased();

	// Output the entire matrix.
	void output(std::ostream &out) const;

	// Iterate through a column.
	template <typename Functor>
	inline void iterX(int x) const;

	// Iterate through a row.
	template <typename Functor>
	inline void iterY(int y) const;

	// Iterate through a column of this matrix and another PreachingBased
	// matrix at the same coordinates.
	template <typename Functor>
	inline void iterX2(int x, PreachingBased &p2) const;

	// Iterate through a row of this matrix and another PreachingBased
	// matrix at the same coordinates.
	template <typename Functor>
	inline void iterY2(int y, PreachingBased &p2) const;

	// Iterate through a column of this matrix and two other PreachingBased
	// matrices at the same coordinates.
	template <typename Functor>
	inline void iterX3(int x, PreachingBased &p2, PreachingBased &p3) const;

	// Iterate through a row of this matrix and two other PreachingBased
	// matrices at the same coordinates.
	template <typename Functor>
	inline void iterY3(int y, PreachingBased &p2, PreachingBased &p3) const;

public:
	// The main data array. This is stored in order of occurrence from left to
	// right, then top to bottom.
	LinkedElm *data;
	// The row pointers. The very last pointer is invalid and is used to denote
	// the end of the array for iterators.
	LinkedElm *Vyc[Z*Y+1];
	// The column pointers.
	LinkedElm *Vxc[Z*X];
	// A reference to the Preaching matrix.
	const Preaching<Y,X,YRHO,XRHO> &preaching;
};


template <typename Elm, int Y, int X, int YRHO, int XRHO>
PreachingBased<Elm,Y,X,YRHO,XRHO>::PreachingBased(const Preaching<Y,X,YRHO,XRHO> &preachingInit) :
	preaching(preachingInit)
{
	// Allocate the memory for the main data array.
	data = new LinkedElm[preaching.ones];

	// Store static pointers to the Vyc and Vxc members for use in the functor.
	static LinkedElm **rVyc, **rVxc;
	rVyc = Vyc;
	rVxc = Vxc;
	// Store a static pointer to the Preaching instance's Hyc array for use in
	// the functor.
	static const int (*rHyc)[Z*Y][YRHO+1];
	rHyc = &preaching.Hyc;

	// A pointer to the current data element.
	LinkedElm *dati = data;
	for (int y = 0;; y++)	// Iterate over the rows.
	{
		rVyc[y] = dati;		// Store the start of the row.
		if (y >= Z*Y)
			break;
		// A pointer to the compressed row in the Preaching matrix.
		const int *pHyc = (*rHyc)[y];
		while (*pHyc++ >= 0)	// Find the end of the row.
			dati++;
	}

	for (int x = 0; x < Z*X; x++)	// Iterate over the columns.
	{
		// Remember the element previous to this one in the column.
		static LinkedElm *prev;
		prev = 0;

		// The functor called for every element in the column.
		struct functor_setlinks {
			static inline void callbackX(int y, int x) {
				// A pointer to the current element in the row.
				LinkedElm *currentelm = rVyc[y];
				// A pointer to the current element in the compressed Preaching
				// matrix row.
				const int *currenth = (*rHyc)[y];
				// Find the element at (y,x).
				while (*currenth != x)
				{
					// Increment both the Preaching and PreachingBased pointers.
					currenth++;
					currentelm++;
				}
				if (prev)						// We're not at the beginning of the column.
					prev->nexty = currentelm;	// Set the previous link.
				else							// We're at the beginning of the column.
					rVxc[x] = currentelm;		// Set the column's beginning.
				prev = currentelm;				// Update the prev pointer.
			}
		};
		// Iterate through every element in the column.
		preaching.iterX<functor_setlinks>(x);
		// Terminate the column.
		prev->nexty = 0;
	}
}

template <typename Elm, int Y, int X, int YRHO, int XRHO>
PreachingBased<Elm,Y,X,YRHO,XRHO>::~PreachingBased()
{
	// Free the main data array.
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
	LinkedElm *e1 = Vyc[y], *e2 = p2.Vyc[y], *e3 = p3.Vyc[y];
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
				out << dati->val << '\t';
				dati++;
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
