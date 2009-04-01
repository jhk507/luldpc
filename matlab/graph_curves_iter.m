% $URL$
% $Date$
% $Rev$

function graph_curves_iter(filename, titlename, axis_snr)
	% Load the histogram data.
	hist = 1 - load(['hist_snr_', filename, '.tsv']);
	
	% Set up the sizes and iteration axis.
	niters = size(hist,2);
	nsiters = 10;
	axis_iter = zeros(nsiters,1);

	% Plot the curves.
	axis_iter_text = cell(nsiters,1);
	for i = 1:nsiters
		axis_iter(i) = max(1, floor(i/nsiters*niters));
		axis_iter_text{i} = ['i=', num2str(axis_iter(i))];

		semilogy(axis_snr, hist(:,axis_iter(i)));
		if i < nsiters
			hold all;
		end
	end
	
	% Create the legend.
	hlegend = legend(axis_iter_text);
	set(hlegend, 'Position', [0 0.5 0.05 0.05]);

	% Make the labels.
	title([titlename, ' block error rate vs. signal-to-noise ratio']);
	xlabel('SNR (dB)');
	ylabel('BLER');
	
	% Set the axes and grids.
	axis tight;
	grid on;
end
