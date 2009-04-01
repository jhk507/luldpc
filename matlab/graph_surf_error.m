% $URL$
% $Date$
% $Rev$

function graph_surf_error(filename, titlename, axis_err, cam, nplot)
	% Load the histogram data.
	hist = load(['hist_err_', filename, '.tsv']);
	axis_iters = 1:size(hist,2);

	% Plot the surface.
	subplot(2,1,nplot);
	surf(axis_iters, axis_err, hist, 'MeshStyle', 'row');
	
	% Set up the axis properties.
	set(gca, 'XScale', 'log');
	set(gca, 'YDir', 'reverse');
	axis tight;

	% Make the labels and colorbar.
	title([titlename, ' error histogram']);
	xlabel('Iteration number');
	ylabel('Error');
	zlabel('Frequency');
	colorbar;
	
	% Set up the camera.
	% campos(cam);
end
