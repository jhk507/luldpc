///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <list>

#include "automatrix.hpp"

#define Z 96	// Matrix expansion factor
#define K (N-M)	// Width of parity portion of Preaching matrix

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Rho is the maximum matrix sparsity parameter

template <unsigned Y, unsigned X, unsigned RHO>
class Preaching
{
public:
	// Constructor initializes compressed matrices
	template <unsigned XH>
	Preaching(const int (&Hinit)[Y][XH], unsigned off);

	// Get a single bit of the expanded matrix
	inline bool at(unsigned y, unsigned x) const;

	// Multiply a row with the expanded matrix

	// Functor member is of the form static void callback(unsigned x, bool p)
	template <typename Functor>
	inline void multRow(const bool (&row)[Z*Y]) const;
	
	inline void multRow(const bool (&row)[Z*Y], bool (&prodrow)[Z*X]) const;

	// Multiply the expanded matrix with a column
	
	// Functor member is of the form inline static void callback(unsigned y, bool p)
	template <typename Functor>
	inline void multCol(const bool (&col)[Z*X]) const;
	
	inline void multCol(const bool (&col)[Z*X], bool (&prodcol)[Z*Y]) const;

	// Iterate over nonzero elements in expanded matrix
	// Both functor members are of the form inline static void callback(unsigned y, unsigned x, unsigned i)
	template <typename Functor>	// For a y, iterate over x
	inline void iterY(unsigned y) const;
	template <typename Functor> // For an x, iterate over y
	inline void iterX(unsigned x) const;

	// Multiply a ZxZ matrix expanded from the Preaching element at (xp,yp),
	// 0 <= xp <= X, 0 <= yp <= Y; with a Zx1 column, and return the index into
	// the OLD column determined by y, the index in the NEW (product) column.
	inline bool pshift(unsigned yp, unsigned xp, const bool *col, unsigned y) const;

public:
	const Automatrix2Alias<int,Y,X> H;	// The matrix array reference

	// These two matrices store the position of "1"s in the expanded matrix,
	// and are terminated with a -1.
	int Hyc[Z*Y][RHO+1];	// Compressed matrix, y-dimension
	int Hxc[Z*X][RHO+1];	// Compressed matrix, x-dimension
	unsigned ones;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


template <unsigned Y, unsigned X, unsigned RHO>
template <unsigned XH>
Preaching<Y,X,RHO>::Preaching(const int (&Hinit)[Y][XH], unsigned off) :
	H(Hinit, off)
{
	ones = 0;

	// Initialize the compressed matrices
	for (unsigned yb = 0; yb < Y; yb++) // Iterate over preaching blocks
	{
		for (unsigned ys = 0; ys < Z; ys++) // Iterate within preaching subblocks
		{
			// Stores the position of "1" elements within expanded Preaching matrix
			int *const pHyc = Hyc[Z*yb+ys];

			unsigned c = 0;
			for (unsigned xb = 0; xb < X; xb++)
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

	for (unsigned xb = 0; xb < X; xb++) // Iterate over preaching blocks
	{
		for (unsigned xs = 0; xs < Z; xs++) // Iterate within preaching subblocks
		{
			// Stores the position of "1" elements within expanded Preaching matrix
			int *const pHxc = Hxc[Z*xb+xs];

			unsigned c = 0;
			for (unsigned yb = 0; yb < Y; yb++)
			{
				const int h = H[yb][xb];
				if (h >= 0)
					pHxc[c++] = (Z - h + xs)%Z + Z*yb;
			}
			pHxc[c] = -1;
		}
	}
}

template <unsigned Y, unsigned X, unsigned RHO>
inline bool Preaching<Y,X,RHO>::at(unsigned y, unsigned x) const
{
/*	// Search for y within Hxc[x]
	const unsigned *pHxc;
	for (pHxc = Hxc[x]; y > *pHxc; pHxc++)
		if (*pHxc < 0)
			return false;
	return *pHxc == y;*/
	const int p = H[y/Z][x/Z];  //
	if (p < 0)
		return false;
	return !((p+y-x)%Z); //converts from unexpanded to expanded bit
}

template <unsigned Y, unsigned X, unsigned RHO>
template <typename Functor>
inline void Preaching<Y,X,RHO>::multRow(const bool (&row)[Z*Y]) const
{
	for (unsigned x = 0; x < Z*X; x++)
	{
		bool prod = 0;
		for (const int *pHxc = Hxc[x]; *pHxc >= 0; pHxc++)
			prod ^= row[*pHxc];
		Functor::callback(x, prod);
	}
}

template <unsigned Y, unsigned X, unsigned RHO>
inline void Preaching<Y,X,RHO>::multRow(const bool (&row)[Z*Y], bool (&prodrow)[Z*X]) const
{
	for (unsigned x = 0; x < Z*X; x++)
	{
		bool prod = 0;
		for (const int *pHxc = Hxc[x]; *pHxc >= 0; pHxc++)
			prod ^= row[*pHxc];
		prodrow[x] = prod;
	}
}

template <unsigned Y, unsigned X, unsigned RHO>
template <typename Functor>
inline void Preaching<Y,X,RHO>::multCol(const bool (&col)[Z*X]) const
{
	for (unsigned y = 0; y < Z*Y; y++)
	{
		bool prod = 0;
		for (const int *pHyc = Hyc[y]; *pHyc >= 0; pHyc++)
			prod ^= col[*pHyc];
		Functor::callback(y, prod);
	}
}

template <unsigned Y, unsigned X, unsigned RHO>
inline void Preaching<Y,X,RHO>::multCol(const bool (&col)[Z*X], bool (&prodcol)[Z*Y]) const
{
	for (unsigned y = 0; y < Z*Y; y++)
	{
		bool prod = 0;
		for (const int *pHyc = Hyc[y]; *pHyc >= 0; pHyc++)
			prod ^= col[*pHyc];
		prodcol[y] = prod;
	}
}

template <unsigned Y, unsigned X, unsigned RHO>
template <typename Functor>
inline void Preaching<Y,X,RHO>::iterY(unsigned y) const
{
	const int *pHyc = Hyc[y];
	int x = *pHyc;
	unsigned i = 0;
	do
	{
		Functor::callbackY(y, i);
		i++;
		pHyc++;
		x = *pHyc;
	} while (x >= 0);
}

template <unsigned Y, unsigned X, unsigned RHO>
template <typename Functor>
inline void Preaching<Y,X,RHO>::iterX(unsigned x) const
{
	const int *pHxc = Hxc[x];
	int y = *pHxc;
	unsigned i = 0;
	do
	{
		Functor::callbackX(x, i);
		i++;
		pHxc++;
		y = *pHxc;
	} while (y >= 0);
}

template <unsigned Y, unsigned X, unsigned RHO>
inline bool Preaching<Y,X,RHO>::pshift(unsigned yp, unsigned xp, const bool *col, unsigned y) const
{
	const int p = H[yp][xp];
	if (p < 0)
		return false;
	return col[(y+p)%Z];
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
