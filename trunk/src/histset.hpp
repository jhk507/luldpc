/*
* $URL$
* $Date$
* $Rev$
*/

#pragma once

#include <fstream>
#include <string>
#include <iomanip>

#include "histogram.hpp"

namespace LDPC
{

template <typename HistType>
class HistogramSet
{
public:
	void init(const char *methodName, const char *errorType,
		const HistType *histsInit)
	{
		hists = histsInit;

		std::string name = methodName;
		name += '_';
		name += errorType;

		std::string filename = "hist_err_";
		filename += name;
		filename += ".tsv";
		ferr.open(filename.c_str());
		ferr << std::setprecision(10);

		filename = "hist_snr_";
		filename += name;
		filename += ".tsv";
		fsnr.open(filename.c_str());
		fsnr << std::setprecision(10);

		filename = "hist_maxiter_";
		filename += name;
		filename += ".tsv";
		fmaxiter.open(filename.c_str());
		fmaxiter << std::setprecision(10);

		filename = "hist_slice_";
		filename += name;
		filename += ".tsv";
		fslice.open(filename.c_str());
		fslice << std::setprecision(10);
	}

	// Output the error histograms.
	void writeLine()
	{
		// Error surface histogram
		// x - iterations
		// y - error buckets
		// z - frequency
		if (snrindex == DEFAULTSNR)
		{
			for (int b = 0; b < NERRBUCKETS; b++)
			{
				for (int i = 0; i < IMAX; i++)
					ferr << hists[i].getNormalizedFreq(b) << '\t';
				ferr << '\n';
			}
		}

		// SNR surface histogram at zero error
		// x - iterations
		// y - SNR
		// z - frequency
		for (int i = 0; i < IMAX; i++)
			fsnr << hists[i].getNormalizedFreq(0) << '\t';
		fsnr << '\n';

		// SNR surface histogram at maximum iteration
		// x - error buckets
		// y - SNR
		// z - frequency
		for (int b = 0; b < NERRBUCKETS; b++)
			fmaxiter << hists[IMAX-1].getNormalizedFreq(b) << '\t';
		fmaxiter << '\n';

		// Giant 4D slice histogram
		// x - iterations
		// y - error buckets
		// y2 - SNR
		// z - frequency
		for (int b = 0; b < NERRBUCKETS; b++)
		{
			for (int i = 0; i < IMAX; i++)
				fslice << hists[i].getNormalizedFreq(b) << '\t';
			fslice << '\n';
		}
	}

private:
	std::ofstream ferr;
	std::ofstream fsnr;
	std::ofstream fmaxiter;
	std::ofstream fslice;

	const HistType *hists;
};

}
