/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include <ostream>

// Output a large column or row matrix of length Z*X; values are tab-separated
template <typename Elm, int X, int Z>
void outputLarge(const Elm *edata, std::ostream &out)
{
	const Elm *pelm = edata;
	for (int x = 0; x < X; x++)
	{
		for (int z = 0; z < Z; z++)
			out << *pelm++ << '\t';
		out << '\t';
	}
	out << std::endl;
}

// Output a large column or row matrix of length Z*X; values are not separated;
// Z blocks are separated.
template <typename Elm, int X, int Z>
void outputLargeContiguous(const Elm *edata, std::ostream &out)
{
	const Elm *pelm = edata;
	for (int x = 0; x < X; x++)
	{
		for (int z = 0; z < Z; z++)
			out << *pelm++;
		out << ' ';
	}
	out << std::endl;
}

// Class to store an array of arbitrary type
template <typename Elm, typename Array>
class Automatrix
{
public:
	// Get a reference to the array
	inline operator Array&()
	{
		return (Array&)*data;
	}

	// Get a pointer to the array data, offset by off
	inline Elm *getData(int off)
	{
		return data+off;
	}

	// Output the array, separated by tabs
	template <int X, int Z>
	inline void outputLarge(std::ostream &out) const
	{
		::outputLarge<Elm,X,Z>(data, out);
	}

	// Output the array; only Z blocks are separated by tabs
	template <int X, int Z>
	inline void outputLargeContiguous(std::ostream &out) const
	{
		::outputLargeContiguous<Elm,X,Z>(data, out);
	}

private:
	// The array data
	Array data;
};

// A single-dimension (row or column) basic matrix
template <typename Elm, int X>
class Automatrix1 : public Automatrix<Elm, Elm[X]>
{
};

// A dual-dimension basic matrix
template <typename Elm, int Y, int X>
class Automatrix2 : public Automatrix<Elm, Elm[Y][X]>
{
};

// An 'alias' into a part of a dual-dimension matrix.
template <typename Elm, int Y, int X>
class Automatrix2Alias
{
public:
	// The type representing a row array.
	typedef const Elm (&SubArrayRef)[X];

public:
	// The constructor. Set the data pointers based on the provided offset.
	template <int XH>
	Automatrix2Alias(const Elm (&original)[Y][XH], int off)
	{
		for (int y = 0; y < Y; y++)
			data[y] = &original[y][off];
	}

	// Make the first [] operator return a reference to a row, suitable for
	// a second [].
	inline SubArrayRef operator[](int y) const
	{
		return (SubArrayRef)*data[y];
	}

private:
	// The array of shifted row pointers.
	const Elm *data[Y];
};
