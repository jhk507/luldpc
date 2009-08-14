#include <fstream>
#include <iomanip>

#include "mathfun.hpp"

using namespace std;

const double tanhlut[] =
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
		flut << tanh((i+0.5)*(TANHMAX/LUTSIZE));
		if (i < LUTSIZE-1)
			flut << ",\n";
		else
			break;
	}
}

