% $URL$
% $Date$
% $Rev$

function graph_curves_snr(filename, titlename, axis_snr)
	% Load the histogram data.
	hist = 1 - load(['hist_snr_', filename, '.tsv']);
	
	% Set up the sizes and iteration axis.
	nsnrs = length(axis_snr);
	niters = size(hist,2);
	axis_iters = 0:niters-1;

	% Plot the curves.
	axis_snr_text = cell(nsnrs,1);
	for s = 1:nsnrs
		axis_snr_text{s} = [num2str(axis_snr(s)),'dB'];

		semilogy(axis_iters, hist(s,:));
		if s < nsnrs
			hold all;
		end
	end
	
	% Create the legend.
	hlegend = legend(axis_snr_text);
	set(hlegend, 'Position', [0 0.5 0.05 0.05]);

	% Make the labels.
	title([titlename, ' block error rate vs. maximum iterations']);
	xlabel('Iteration number');
	ylabel('BLER');
	
	% Set the axes and grids.
	axis tight;
	grid on;
end
