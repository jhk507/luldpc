/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include <ostream>

// Output a large column or row matrix of length Z*X; values are tab-separated
template <int X, int Z, typename Elm>
void outputLarge(const Elm (&edata)[Z*X], std::ostream &out)
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
template <int X, int Z, typename Elm>
void outputLargeContiguous(const Elm (&edata)[Z*X], std::ostream &out)
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
