% $URL$
% $Date$
% $Rev$

function graph_curves_perf(axistitle)
	% Declare the global variables.
	global ndecodes axis_decode;
	
	% Load the histogram data.
	hist = load('hist_perf.tsv');
	
	% Load the performance axis.
	axis_perf = load('axis_perf.tsv');
	
	% Plot the curves.
	for d = 1:ndecodes
		plot(axis_perf, hist(d,:));
		if d < ndecodes
			hold all;
		end
	end
	
	% Make the labels and legend.
	title(axistitle);
	xlabel('Iteration duration (us)');
	ylabel('Frequency');
	set(legend(axis_decode, 'Location', 'BestOutside'), 'Interpreter','none');
	
	% Set the axes and grids.
	axis tight;
	grid on;
end
