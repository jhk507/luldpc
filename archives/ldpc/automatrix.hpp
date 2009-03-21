#ifndef INCLUDED_AUTOMATRIX_HPP
#define INCLUDED_AUTOMATRIX_HPP


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

template <typename Elm, unsigned X>
class Automatrix1 : public Automatrix<Elm, Elm[X]>
{
};

template <typename Elm, unsigned Y, unsigned X>
class Automatrix2 : public Automatrix<Elm, Elm[Y][X]>
{
};

template <typename Elm, unsigned Y, unsigned X>
class Automatrix2Alias
{
public:
	typedef const Elm (&SubArrayRef)[X];

public:
	template <unsigned XH>
	Automatrix2Alias(const Elm (&original)[Y][XH], unsigned off);

	inline SubArrayRef operator[](unsigned y) const
	{
		return (SubArrayRef)*data[y];
	}

private:
	const Elm *data[Y];
};

template <typename Elm, unsigned Y, unsigned X>
template <unsigned XH>
Automatrix2Alias<Elm,Y,X>::Automatrix2Alias(const Elm (&original)[Y][XH], unsigned off)
{
	for (unsigned y = 0; y < Y; y++)
		data[y] = &original[y][off];
}


#endif
