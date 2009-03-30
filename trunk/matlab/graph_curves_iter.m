% $URL$
% $Date$
% $Rev$

function graph_curves_iter(filename, titlename, axis_snr, nplot)
	% Load the histogram data.
	hist = 1 - load(['hist_snr_', filename, '.tsv']);
	
	% Set up the sizes and iteration axis.
	niters = size(hist,2);
	axis_iter = zeros(10,1);

	% Plot the curves.
	subplot(2,1,nplot);
	for i = 1:length(axis_iter)
		axis_iter(i) = max(1, ...
			floor(i/length(axis_iter)*niters));
		axis_iter_text(i,:) = ['i=', num2str(axis_iter(i),'%0.3d')];

		semilogy(axis_snr, hist(:,axis_iter(i)));
		if (i < length(axis_iter))
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
	grid on;
	axis tight;
end
