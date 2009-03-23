/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include "automatrix.hpp"

#define Z 96	// Matrix expansion factor
#define K (N-M)	// Width of parity portion of Preaching matrix

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Rho is the maximum matrix sparsity parameter

template <int Y, int X, int YRHO, int XRHO>
class Preaching
{
public:
	// Constructor initializes compressed matrices
	template <int XH>
	Preaching(const int (&Hinit)[Y][XH], int off);

	// Get a single bit of the expanded matrix
	inline bool at(int y, int x) const;

	// Multiply a row with the expanded matrix

	// Functor member is of the form static void callback(int x, bool p)
	template <typename Functor>
	inline void multRow(const bool (&row)[Z*Y]) const;
	
	inline void multRow(const bool (&row)[Z*Y], bool (&prodrow)[Z*X]) const;

	// Multiply the expanded matrix with a column
	
	// Functor member is of the form inline static void callback(int y, bool p)
	template <typename Functor>
	inline void multCol(const bool (&col)[Z*X]) const;
	
	inline void multCol(const bool (&col)[Z*X], bool (&prodcol)[Z*Y]) const;

	// Iterate over nonzero elements in expanded matrix
	
	// Functor member is of the form inline static void callbackY(int y, int x, int xi)
	template <typename Functor>	// For a y, iterate over x
	inline void iterY(int y) const;
	// Functor member is of the form inline static void callbackX(int y, int x, int yi)
	template <typename Functor> // For an x, iterate over y
	inline void iterX(int x) const;

	// Multiply a ZxZ matrix expanded from the Preaching element at (xp,yp),
	// 0 <= xp <= X, 0 <= yp <= Y; with a Zx1 column, and return the index into
	// the OLD column determined by y, the index in the NEW (product) column.
	inline bool pshift(int yp, int xp, const bool *col, int y) const;

	void output(std::ostream &out) const;

public:
	const Automatrix2Alias<int,Y,X> H;	// The matrix array reference

	// These two matrices store the position of "1"s in the expanded matrix,
	// and are terminated with a -1.
	int Hyc[Z*Y][YRHO+1];	// Compressed matrix, y-dimension
	int Hxc[Z*X][XRHO+1];	// Compressed matrix, x-dimension
	int ones;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


template <int Y, int X, int YRHO, int XRHO>
template <int XH>
Preaching<Y,X,YRHO,XRHO>::Preaching(const int (&Hinit)[Y][XH], int off) :
	H(Hinit, off)
{
	ones = 0;

	// Initialize the compressed matrices
	for (int yb = 0; yb < Y; yb++) // Iterate over preaching blocks
	{
		for (int ys = 0; ys < Z; ys++) // Iterate within preaching subblocks
		{
			// Stores the position of "1" elements within expanded Preaching matrix
			int *const pHyc = Hyc[Z*yb+ys];

			int c = 0;
			for (int xb = 0; xb < X; xb++)
			{
				const int h = H[yb][xb];
				if (h >= 0)
				{
					pHyc[c++] = (h + ys)%Z + Z*xb;
					ones++;
				}
			}
			pHyc[c] = -1;
		}
	}

	for (int xb = 0; xb < X; xb++) // Iterate over preaching blocks
	{
		for (int xs = 0; xs < Z; xs++) // Iterate within preaching subblocks
		{
			// Stores the position of "1" elements within expanded Preaching matrix
			int *const pHxc = Hxc[Z*xb+xs];

			int c = 0;
			for (int yb = 0; yb < Y; yb++)
			{
				const int h = H[yb][xb];
				if (h >= 0)
					pHxc[c++] = (Z - h + xs)%Z + Z*yb;
			}
			pHxc[c] = -1;
		}
	}
}

template <int Y, int X, int YRHO, int XRHO>
inline bool Preaching<Y,X,YRHO,XRHO>::at(int y, int x) const
{
/*	// Search for y within Hxc[x]
	const int *pHxc;
	for (pHxc = Hxc[x]; y > *pHxc; pHxc++)
		if (*pHxc < 0)
			return false;
	return *pHxc == y;*/
	const int p = H[y/Z][x/Z];  //
	if (p < 0)
		return false;
	return !((p+y-x)%Z); //converts from unexpanded to expanded bit
}

template <int Y, int X, int YRHO, int XRHO>
template <typename Functor>
inline void Preaching<Y,X,YRHO,XRHO>::multRow(const bool (&row)[Z*Y]) const
{
	for (int x = 0; x < Z*X; x++)
	{
		bool prod = 0;
		for (const int *pHxc = Hxc[x]; *pHxc >= 0; pHxc++)
			prod ^= row[*pHxc];
		Functor::callback(x, prod);
	}
}

template <int Y, int X, int YRHO, int XRHO>
inline void Preaching<Y,X,YRHO,XRHO>::multRow(const bool (&row)[Z*Y], bool (&prodrow)[Z*X]) const
{
	for (int x = 0; x < Z*X; x++)
	{
		bool prod = 0;
		for (const int *pHxc = Hxc[x]; *pHxc >= 0; pHxc++)
			prod ^= row[*pHxc];
		prodrow[x] = prod;
	}
}

template <int Y, int X, int YRHO, int XRHO>
template <typename Functor>
inline void Preaching<Y,X,YRHO,XRHO>::multCol(const bool (&col)[Z*X]) const
{
	for (int y = 0; y < Z*Y; y++)
	{
		bool prod = 0;
		for (const int *pHyc = Hyc[y]; *pHyc >= 0; pHyc++)
			prod ^= col[*pHyc];
		Functor::callback(y, prod);
	}
}

template <int Y, int X, int YRHO, int XRHO>
inline void Preaching<Y,X,YRHO,XRHO>::multCol(const bool (&col)[Z*X], bool (&prodcol)[Z*Y]) const
{
	for (int y = 0; y < Z*Y; y++)
	{
		bool prod = 0;
		for (const int *pHyc = Hyc[y]; *pHyc >= 0; pHyc++)
			prod ^= col[*pHyc];
		prodcol[y] = prod;
	}
}

template <int Y, int X, int YRHO, int XRHO>
template <typename Functor>
inline void Preaching<Y,X,YRHO,XRHO>::iterY(int y) const
{
	const int *pHyc = Hyc[y];
	int x = *pHyc;
	int xi = 0;
	do
	{
		Functor::callbackY(y, x, xi);
		xi++;
		pHyc++;
		x = *pHyc;
	} while (x >= 0);
}

template <int Y, int X, int YRHO, int XRHO>
template <typename Functor>
inline void Preaching<Y,X,YRHO,XRHO>::iterX(int x) const
{
	const int *pHxc = Hxc[x];
	int y = *pHxc;
	int yi = 0;
	do
	{
		Functor::callbackX(y, x, yi);
		yi++;
		pHxc++;
		y = *pHxc;
	} while (y >= 0);
}

template <int Y, int X, int YRHO, int XRHO>
inline bool Preaching<Y,X,YRHO,XRHO>::pshift(int yp, int xp, const bool *col, int y) const
{
	const int p = H[yp][xp];
	if (p < 0)
		return false;
	return col[(y+p)%Z];
}

template <int Y, int X, int YRHO, int XRHO>
void Preaching<Y,X,YRHO,XRHO>::output(std::ostream &out) const
{
	for (int y = 0; y < Y*Z; y++)
	{
		for (int x = 0; x < X*Z; x++)
		{
			out << at(y,x);
			if (x%Z == Z-1)
				out << ' ';
		}
		out << endl;
		if (y%Z == Z-1)
			out << endl;
	}
	out << endl;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
