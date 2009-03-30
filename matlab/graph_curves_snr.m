% $URL$
% $Date$
% $Rev$

function graph_curves_snr(filename, titlename, axis_snr, nplot)
	% Load the histogram data.
	hist = load(['hist_snr_', filename, '.tsv']);
	
	% Set up the sizes and iteration axis.
	nsnrs = size(hist,1);
	niters = size(hist,2);
	axis_iters = 0:niters-1;
	
	% Plot the curves.
	subplot(2,1,nplot);
	for s = 1:nsnrs
		semilogy(axis_iters, 1 - hist(s,:));
		hold on;
	end

	% Make the labels.
	title([titlename, ' block error rate vs. maximum iterations']);
	xlabel('Iteration number');
	ylabel('BLER');
	
	grid on;
	axis tight;
end
