% $URL$
% $Date$
% $Rev$

function graph_surf_snr(filename, axistitle, nplot)
	% Load the histogram data.
	hist = load(['hist_snr_', filename, '.tsv']);

	% Declare the global variables.
	global axis_iter;
	global axis_snr;
	
	% Plot the surface.
	subplot(1,2,nplot);
	surf(axis_iter, axis_snr, hist, 'MeshStyle', 'row');
	
	% Set up the axis properties.
	set(gca, 'XScale', 'log');
	axis tight;

	% Make the labels.
	title(axistitle);
	xlabel('Iteration number');
	ylabel('SNR (dB)');
	zlabel('Frequency');
end
