#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

namespace LDPC
{

const double tanhlut[1024] =
{
#include "tanhlut.cpp"
};

void makeluts()
{
	ofstream flut;
	flut << setprecision(30);

	flut.open("tanhlut.cpp");

	for (int i = 0;; i++)
	{
		flut << tanh((i+0.5)*(19.1/1024.0));
		if (i < 1023)
			flut << ",\n";
		else
			break;
	}
}

double tanhapp(double x)
{
	const double xabs = fabs(x);
	const double pval = (xabs >= 19.1) ? 1 :
		tanhlut[(int)(xabs*(1024.0/19.1))];
	return (x >= 0) ? pval : -pval;
}

}
