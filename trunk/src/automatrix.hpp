#pragma once

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
