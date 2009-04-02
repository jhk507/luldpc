% $URL$
% $Date$
% $Rev$

function graph_curves_iter(filename, axistitle)
	% Load the histogram data.
	hist = 1 - load(['hist_snr_', filename, '.tsv']);
	
	% Declare the global variables.
	global nsiters;
	global axis_siter;
	global axis_siter_text;
	global axis_snr;
	
	% Plot the curves.
	for i = 1:nsiters
		semilogy(axis_snr, hist(:,axis_siter(i)));
		if i < nsiters
			hold all;
		end
	end
	
	% Make the labels and legend.
	title(axistitle);
	xlabel('SNR (dB)');
	ylabel('BLER');
	legend(axis_siter_text, 'Location', 'BestOutside');
	
	% Set the axes and grids.
	axis tight;
	grid on;
end
