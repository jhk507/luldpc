% $URL$
% $Date$
% $Rev$

function graph_curves_aveiter(filename, axistitle)
	% Load the histogram data.
	hist = load(['hist_snr_', filename, '.tsv']);
	
	% Declare the global variables.
	global niters;
	global nsnrs;
	global axis_iter;
	global axis_snr_text;
	
	% Convert from a cumulative histogram to a moving weighted average
	wsum = zeros(nsnrs, 1);
	aves = zeros(size(hist));
	for i = 1:niters
		weight = hist(:,i);
		if i > 1
			weight = weight - hist(:,i-1);
		end
		weight = weight*i;
		wsum = wsum + weight;
		weight = wsum;
		for s = 1:nsnrs
			if hist(s,i) ~= 0
				weight(s) = weight(s)./hist(s,i);
			end
		end
		aves(:,i) = weight;
	end

	% Plot the curves.
	for s = 1:nsnrs
		semilogx(axis_iter, aves(s,:));
		if s < nsnrs
			hold all;
		end
	end
	
	% Make the labels and legend.
	title(axistitle);
	xlabel('Maximum iterations');
	ylabel('Average iterations to resolution');
	legend(axis_snr_text, 'Location', 'BestOutside');
	
	% Set the axes and grids.
	axis tight;
	grid on;
end
