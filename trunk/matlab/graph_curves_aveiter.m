% $URL$
% $Date$
% $Rev$

function graph_curves_aveiter(filename, axistitle, axis_snr)
	% Load the histogram data.
	hist = load(['hist_snr_', filename, '.tsv']);
	
	% Set up the sizes and iteration axis.
	nsnrs = length(axis_snr);
	niters = size(hist,2);
	axis_iters = 1:niters;
	
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
	axis_snr_text = cell(nsnrs, 1);
	for s = 1:nsnrs
		axis_snr_text{s} = [num2str(axis_snr(s)),'dB'];

		semilogx(axis_iters, aves(s,:));
		if s < nsnrs
			hold all;
		end
	end
	
	% Create the legend.
	hlegend = legend(axis_snr_text);
	set(hlegend, 'Position', [0.01 0.75 0.01 0.01]);

	% Make the labels.
	title(axistitle);
	xlabel('Maximum iterations');
	ylabel('Average iterations to resolution');
	
	% Set the axes and grids.
	axis tight;
	grid on;
end
