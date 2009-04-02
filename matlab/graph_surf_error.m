% $URL$
% $Date$
% $Rev$

function graph_surf_error(filename, axistitle, axis_err, cam, nplot)
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

	% Make the labels.
	title(axistitle);
	xlabel('Iteration number');
	ylabel('Error');
	zlabel('Frequency');
	
	% Set up the camera.
	% campos(cam);
end
