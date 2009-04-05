% $URL$
% $Date$
% $Rev$

function graph_surf_maxiter(filename, axistitle, axis_err, cam, nplot)
	% Load the histogram data.
	hist = load(['hist_maxiter_', filename, '.tsv']);
	
	% Declare the global variables.
	global axis_snr;

	% Plot the surface.
	subplot(1,2,nplot);
	surf(axis_err, axis_snr, hist, 'MeshStyle', 'row');
	
	% Set up the axis properties.
	axis tight;

	% Make the labels.
	title(axistitle);
	xlabel('Error');
	ylabel('SNR (dB)');
	zlabel('Frequency');

	% Set up the camera.
	campos(cam);
end
