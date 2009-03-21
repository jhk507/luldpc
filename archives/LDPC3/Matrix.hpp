#pragma once;
/*
template <typename Type, unsigned xm, unsigned ym>
class Matrix
{
protected:
	Matrix();
public:
	virtual ~Matrix();

	virtual Type get(unsigned x, unsigned y) const = 0;
	virtual void set(unsigned x, unsigned y, typename Type) const = 0;

	template <unsigned xm1>
	virtual Matrix<Type, xm1, ym> *operator*(Matrix<Type, xm1, xm> &m1) const
	{
		Matrix<Type, xm1, ym> *const p = new Matrix<Type, xm1, ym>();

		for (unsigned y0 = 0; y0 < ym; y0++)
		{
			for (unsigned x1 = 0; x1 < xm1; x1++)
			{
				Type sum = 0;
				for (unsigned x0 = 0; x0 < xm; x0++)
					sum += get(x0, y0) * m1.get(x1, x0);
				p->set(x1, y0);
			}
		}
	}
};
*/