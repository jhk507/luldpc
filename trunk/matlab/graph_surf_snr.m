% $URL$
% $Date$
% $Rev$

function graph_surf_snr(filename, axistitle, axis_snr, nplot)
	% Load the histogram data.
	hist = load(['hist_snr_', filename, '.tsv']);
	axis_iters = 1:size(hist,2);

	% Plot the surface.
	subplot(2,1,nplot);
	surf(axis_iters, axis_snr, hist, 'MeshStyle', 'row');
	
	% Set up the axis properties.
	set(gca, 'XScale', 'log');
	axis tight;

	% Make the labels.
	title(axistitle);
	xlabel('Iteration number');
	ylabel('SNR (dB)');
	zlabel('Frequency');
end
