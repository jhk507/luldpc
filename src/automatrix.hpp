#pragma once

#include <ostream>


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


template <typename Elm, typename Array>
class Automatrix
{
protected:
	Automatrix()
	{
	}

public:
	inline operator Array&()
	{
		return (Array&)*data;
	}

	inline Elm *getData(int off)
	{
		return data+off;
	}

	template <int X, int Z>
	void outputLarge(std::ostream &out) const
	{
		::outputLarge<Elm,X,Z>(data, out);
	}

	template <int X, int Z>
	void outputLargeContiguous(std::ostream &out) const
	{
		::outputLargeContiguous<Elm,X,Z>(data, out);
	}



protected:
	Array data;
};

template <typename Elm, int X>
class Automatrix1 : public Automatrix<Elm, Elm[X]>
{
};

template <typename Elm, int Y, int X>
class Automatrix2 : public Automatrix<Elm, Elm[Y][X]>
{
};

template <typename Elm, int Y, int X>
class Automatrix2Alias
{
public:
	typedef const Elm (&SubArrayRef)[X];

public:
	template <int XH>
	Automatrix2Alias(const Elm (&original)[Y][XH], int off)
	{
		for (int y = 0; y < Y; y++)
			data[y] = &original[y][off];
	}


	inline SubArrayRef operator[](int y) const
	{
		return (SubArrayRef)*data[y];
	}

private:
	const Elm *data[Y];
};
