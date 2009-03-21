#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>

#include "MTRand_gaussian.hpp"
#include "PreachingFast.hpp"
#include "PreachingReference.hpp"

using namespace std;


int main()
{
	const int XM = 24, YM = 12;

	// Preaching matrix H
	int H[YM][XM] =
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

	// Discrete value random number generator
	const unsigned long seed = (unsigned long)time(NULL);
	MTRand_int32 rand(1224098149);
	// Random bits
	bool *const r = new char[XM*96];
	// Intermediary sum of random bits
	bool rs[XM] = { 0 };
	// Generate the random bits and set the intermediary sum
	for (int i = 0; i < XM*96; i++)
	{
		const bool rn = rand() & 1;
		r[i] = rn;
		rs[i%96] ^= rn;
	}






	/************************************************************ // Testing the encoder
	*************************************************************/
	bool *const m12 = new bool[1152*2304];

	for (int Y = 0; Y < 12; Y++)
	{
		for (int X = 0; X < 24; X++)
		{
			const int M = H[Y][X];
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

	bool *const cresult = new bool[1152];
	for (int y = 0; y < 1152; y++)
	{
		cresult[y] = 0;
		for (int x = 0; x < 2304; x++)
			cresult[y] ^= r[x] & m12[y*2304+x];
	}




	// Do the matrix multiplication of the Preaching matrix and the random matrix
	// Uses sparse shortcuts.
	bool *const x = new bool[YM*96];
	for (int Y = 0; Y < YM; Y++)
	{
		for (int y = 0; y < 96; y++)
		{
			bool xresult = 0;
			for (int X = 0; X < XM; X++)
			{
				// The current element in the preaching matrix.
				const int h = H[Y][X];
				switch (h)
				{
				case -1: // This row is all 0, so there's nothing to add
					break;
				case 0: // This row is all 1, so add the precalculated sum
					xresult ^= rs[X];
					break;
				default: // This row only has two '1' values at predefined positions
					if (y < h)
					{
						const bool *const rp = r + X*96 + h-y-1;
						xresult ^= (*rp) ^ rp[1];
					}
					else if (y == h)
						xresult ^= r[X*96];
				}
			}
			x[Y*96+y] = xresult;
		}
	}


	bool equal = true;
	for (int i = 0; i < 1152; i++)
	{
		if (x[i] != cresult[i])
		{
			equal = false;
			break;
		}
	}

	delete [] m12;
	

	double *const y = new double[1152]; // The signal
	// The Gaussian distribution random number generator
	MTRand_gaussian grand((unsigned long)time(NULL));

	// Constants for AWGN calculation
	const double R = 0.5;
	const double snr = pow(10.0,1.5/10);
	const double sigma = pow(2*R*snr, -0.5);
	
	// Perform BPSK and AWGN addition
	for (int i = 0; i < XM*96; i++)
		y[i] = (x[i] ? 1 : -1) + grand()*sigma; //use ternary operator

	delete [] r;
	delete [] cresult;
	delete [] y;
	return 0;
}

/*	const int pts[3][2] = { //TESTING SPARSE MATRIX
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

