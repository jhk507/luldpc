#include <time.h>
#include <iostream>
#include <fstream>

#include "mtrand/mtrand.hpp"

using namespace std;

int main()
{
	int m12p[12][24] =
	{
		{-1,94,73,-1,-1,-1,-1,-1,55,83,-1,-1, 7, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
		{-1,27,-1,-1,-1,22,79, 9,-1,-1,-1,12,-1, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1},
		{-1,-1,-1,24,22,81,-1,33,-1,-1,-1, 0,-1,-1, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1},
		{61,-1,47,-1,-1,-1,-1,-1,65,25,-1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,-1,-1},
		{-1,-1,39,-1,-1,-1,84,-1,-1,41,72,-1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,-1},
		{-1,-1,-1,-1,46,40,-1,82,-1,-1,-1,79, 0,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1},
		{-1,-1,95,53,-1,-1,-1,-1,-1,14,18,-1,-1,-1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1},
		{-1,11,73,-1,-1,-1, 2,-1,-1,47,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,-1,-1,-1},
		{12,-1,-1,-1,83,24,-1,43,-1,-1,-1,51,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,-1,-1},
		{-1,-1,-1,-1,-1,94,-1,59,-1,-1,70,72,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,-1},
		{-1,-1, 7,65,-1,-1,-1,-1,39,49,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0},
		{43,-1,-1,-1,-1,66,-1,41,-1,-1,-1,26, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0}
	};

	MTRand_int32 rand(time(NULL));
	bool rbits[2304];
	bool rbitsp[24] = { 0 };
	for (int i = 0; i < 2304; i++)
	{
		const bool r = rand() & 1;
		rbits[i] = r;
		rbitsp[i%96] ^= r;
	}

	bool result[1152];
	for (int Y = 0; Y < 12; Y++)
	{
		for (int y = 0; y < 96; y++)
		{
			bool nresult = 0;
			for (int X = 0; X < 24; X++)
			{
				const int M = m12p[Y][X];
				switch (M)
				{
				case -1:
					break;
				case 0:
					nresult ^= rbitsp[X];
					break;
				default:
					if (y < M)
					{
						const bool *const r = rbits + X*96 + M-y-1;
						nresult ^= r[0] ^ r[1];
					}
					else if (y == M)
						nresult ^= rbits[X*96];
				}
			}
			result[Y*96+y] = nresult;
		}
	}

	/************************************************************
	*************************************************************
	bool *m12 = new bool[1152*2304];

	for (int Y = 0; Y < 12; Y++)
	{
		for (int X = 0; X < 24; X++)
		{
			const int M = m12p[Y][X];
			switch (M)
			{
				case 0:
				for (int y = 0; y < 96; y++)
					for (int x = 0; x < 96; x++)
						m12[(Y*96+y)*2304+X*96+x] = 1;
				break;
				case -1:
				for (int y = 0; y < 96; y++)
					for (int x = 0; x < 96; x++)
						m12[(Y*96+y)*2304+X*96+x] = 0;
				break;
				default:
				for (int y = 0; y < 96; y++)
					for (int x = 0; x < 96; x++)
						m12[(Y*96+y)*2304+X*96+x] = (x == M - y) || (x == M - y - 1);
			}
		}
	}

	bool cresult[1152];
	for (int y = 0; y < 1152; y++)
	{
		cresult[y] = 0;
		for (int x = 0; x < 2304; x++)
			cresult[y] ^= rbits[x] & m12[y*2304+x];
	}
	bool equal = true;
	for (int i = 0; i < 1152; i++)
	{
		if (result[i] != cresult[i])
		{
			equal = false;
			break;
		}
	}

	delete [] m12;
	*/

	return 0;
}

/*	const int pts[3][2] = {
		{0,0},
		{0,12},
		{0,13}
	};

	ofstream testo("testoutput.txt");


	for (int n = 0; n < 3; n++)
	{
		const int X = pts[n][1];
		const int Y = pts[n][0];
		testo << '(' << X << ',' << Y << "): " << m12p[Y][X] << endl;
		for (int y = 0; y < 96; y++)
		{
			for (int x = 0; x < 96; x++)
				testo << (unsigned)m12[(Y*96+y)*2304+X*96+x];
			testo << endl;
		}
		testo << endl;
	}

	testo << "Result of random matrix multiplication: " << endl;
	for (int n = 0; n < 1152; n++)
		testo << result[n];
	testo << endl;

	delete [] m12;
*/

