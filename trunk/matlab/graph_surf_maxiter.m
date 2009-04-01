% $URL$
% $Date$
% $Rev$

function graph_surf_maxiter(filename, titlename, axis_err, axis_snr, cam, nplot)
	% Load the histogram data.
	hist = load(['hist_maxiter_', filename, '.tsv']);

	% Plot the surface.
	subplot(2,1,nplot);
	surf(axis_err, axis_snr, hist, 'MeshStyle', 'row');
	
	% Set up the axis properties.
	axis tight;

	% Make the labels and colorbar.
	title([titlename, ' maximal iteration error histogram']);
	xlabel('Error');
	ylabel('SNR (dB)');
	zlabel('Frequency');
	colorbar;

	% Set up the camera.
	campos(cam);
end
