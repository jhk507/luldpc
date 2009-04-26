/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include "automatrix.hpp"

#define Z 96				// Matrix expansion factor
#define K (N-M)				// Width of parity portion of Preaching matrix
#define RATE ((double)M/N)	// Rate

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// A Preaching class that constructs a sparse matrix
// Rho is the maximum matrix sparsity parameter for the X and Y dimensions.
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

	// Functor member is of the form static void callbackProduct(int x, bool p)
	template <typename Functor>
	inline void multRowCallback(Functor &func, const bool (&row)[Z*Y]) const;

	inline void multRow(const bool (&row)[Z*Y], bool (&prodrow)[Z*X]) const;

	// Multiply the expanded matrix with a column

	// Functor member is of the form inline static void callbackProduct(int y, bool p)
	template <typename Functor>
	inline void multColCallback(Functor &func, const bool (&col)[Z*X]) const;

	inline void multCol(const bool (&col)[Z*X], bool (&prodcol)[Z*Y]) const;

	// Iterate over nonzero elements in expanded matrix

	// Functor member is of the form inline static void callbackY(int y, int x)
	template <typename Functor>	// For a y, iterate over x
	inline void iterY(Functor &func, int y) const;
	// Functor member is of the form inline static void callbackX(int y, int x)
	template <typename Functor> // For an x, iterate over y
	inline void iterX(Functor &func, int x) const;

	// Multiply a ZxZ matrix expanded from the Preaching element at (xp,yp),
	// 0 <= xp <= X, 0 <= yp <= Y; with a Zx1 column, and return the index into
	// the OLD column determined by y, the index in the NEW (product) column.
	inline bool pshift(int yp, int xp, const bool *col, int y) const;

	// Output the expanded Preaching matrix.
	void output(std::ostream &out) const;

public:
	// The unexpanded matrix array reference
	const Automatrix2Alias<int,Y,X> H;

	// These two matrices store the position of "1"s in the expanded matrix,
	// and are terminated with a -1.
	int Hyc[Z*Y][YRHO+1];	// Compressed matrix, y-dimension
	int Hxc[Z*X][XRHO+1];	// Compressed matrix, x-dimension
	// The number of ones in the expanded matrix.
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
	for (int yb = 0; yb < Y; yb++) // Iterate over preaching blocks, y
	{
		for (int ys = 0; ys < Z; ys++) // Iterate within preaching subblocks, y
		{
			// Set a pointer to the current compressed y-dimension data
			int *const pHyc = Hyc[Z*yb+ys];

			int c = 0;
			for (int xb = 0; xb < X; xb++) // Iterate over preaching blocks, x
			{
				// Examine the current unexpanded value
				const int h = H[yb][xb];
				if (h >= 0)
				{
					// Stores the position of the "1" element within expanded Preaching matrix
					pHyc[c++] = (h + ys)%Z + Z*xb;
					ones++;
				}
			}
			// Terminate the compressed data
			pHyc[c] = -1;
		}
	}

	for (int xb = 0; xb < X; xb++)			// Iterate over preaching blocks, x
	{
		for (int xs = 0; xs < Z; xs++)		// Iterate within preaching subblocks, x
		{
			// Set a pointer to the current compressed x-dimension data
			int *const pHxc = Hxc[Z*xb+xs];

			int c = 0;
			for (int yb = 0; yb < Y; yb++)	// Iterate over preaching blocks, y
			{
				// Examine the current unexpanded value
				const int h = H[yb][xb];
				if (h >= 0)
					// Stores the position of the "1" element within expanded Preaching matrix
					pHxc[c++] = (Z - h + xs)%Z + Z*yb;
			}
			// Terminate the compressed data
			pHxc[c] = -1;
		}
	}
}

template <int Y, int X, int YRHO, int XRHO>
inline bool Preaching<Y,X,YRHO,XRHO>::at(int y, int x) const
{
	const int p = H[y/Z][x/Z];  // Get the appropriate unexpanded value
	if (p < 0)					// If negative, the whole block is 0
		return false;
	return !((p+y-x)%Z);		// Checks to see if there's a 1 at these coordinates
}

template <int Y, int X, int YRHO, int XRHO>
template <typename Functor>
inline void Preaching<Y,X,YRHO,XRHO>::multRowCallback(Functor &func, const bool (&row)[Z*Y]) const
{
	// Do the matrix multiplication loop, over x first
	for (int x = 0; x < Z*X; x++)
	{
		bool prod = 0;
		for (const int *pHxc = Hxc[x]; *pHxc >= 0; pHxc++)
			prod ^= row[*pHxc];
		func.callbackProduct(x, prod);
	}
}

template <int Y, int X, int YRHO, int XRHO>
inline void Preaching<Y,X,YRHO,XRHO>::multRow(const bool (&row)[Z*Y], bool (&prodrow)[Z*X]) const
{
	// Do the matrix multiplication loop, over x first
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
inline void Preaching<Y,X,YRHO,XRHO>::multColCallback(Functor &func, const bool (&col)[Z*X]) const
{
	// Do the matrix multiplication loop, over y first
	for (int y = 0; y < Z*Y; y++)
	{
		bool prod = 0;
		for (const int *pHyc = Hyc[y]; *pHyc >= 0; pHyc++)
			prod ^= col[*pHyc];
		func.callbackProduct(y, prod);
	}
}

template <int Y, int X, int YRHO, int XRHO>
inline void Preaching<Y,X,YRHO,XRHO>::multCol(const bool (&col)[Z*X], bool (&prodcol)[Z*Y]) const
{
	// Do the matrix multiplication loop, over y first
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
inline void Preaching<Y,X,YRHO,XRHO>::iterY(Functor &func, int y) const
{
	const int *pHyc = Hyc[y];	// A pointer to the current row.
	int x = *pHyc;				// The current x value.
	// Loop through the row, calling the callback for every coordinate present
	// in the sparse matrix.
	do
	{
		func.callbackY(y, x);
		pHyc++;
		x = *pHyc;
	} while (x >= 0);
}

template <int Y, int X, int YRHO, int XRHO>
template <typename Functor>
inline void Preaching<Y,X,YRHO,XRHO>::iterX(Functor &func, int x) const
{
	const int *pHxc = Hxc[x];	// A pointer to the current column.
	int y = *pHxc;				// The current y value.
	do
	{
		// Loop through the column, calling the callback for every coordinate
		// present in the sparse matrix.
		func.callbackX(y, x);
		pHxc++;
		y = *pHxc;
	} while (y >= 0);
}

template <int Y, int X, int YRHO, int XRHO>
inline bool Preaching<Y,X,YRHO,XRHO>::pshift(int yp, int xp, const bool *col, int y) const
{
	const int p = H[yp][xp];	// The associated unexpanded value.
	if (p < 0)					// If it's -1, the whole block is zero.
		return false;
	return col[(y+p)%Z];		// Calculate the shifted value.
}

template <int Y, int X, int YRHO, int XRHO>
void Preaching<Y,X,YRHO,XRHO>::output(std::ostream &out) const
{
	// Loop through the entire expanded matrix, outputting a '1' when appropriate.
	for (int y = 0; y < Y*Z; y++)
	{
		for (int x = 0; x < X*Z; x++)
		{
			out << at(y,x);
			if (x%Z == Z-1)
				out << ' ';
		}
		out << std::endl;
		if (y%Z == Z-1)
			out << std::endl;
	}
	out << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
