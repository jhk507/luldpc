% $URL$
% $Date$
% $Rev$

function graph_curves_perf(axistitle)
	% Load the histogram data.
	hist = load('hist_perf.tsv');
	
	% Load the performance axis.
	axis_perf = load('axis_perf.tsv');
	
	% Declare the global variables.
	global ndecodes;
	global axis_decode;
	
	% Plot the curves.
	for d = 1:ndecodes
		plot(axis_perf, hist(d,:));
		if d < ndecodes
			hold all;
		end
	end
	
	% Make the labels and legend.
	title(axistitle);
	xlabel('Average iteration duration (us)');
	ylabel('Frequency');
	legend(axis_decode, 'Location', 'BestOutside');
	
	% Set the axes and grids.
	axis tight;
	grid on;
end
