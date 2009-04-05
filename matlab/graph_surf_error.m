% $URL$
% $Date$
% $Rev$

function graph_surf_error(filename, axistitle, axis_err, nplot)
	% Load the histogram data.
	hist = load(['hist_err_', filename, '.tsv']);
	
	% Declare the global variables.
	global axis_iter;

	% Plot the surface.
	subplot(1,2,nplot);
	surf(axis_iter, axis_err, hist, 'MeshStyle', 'row');
	
	% Set up the axis properties.
	set(gca, 'XScale', 'log', 'YDir', 'reverse');
	axis tight;

	% Make the labels.
	title(axistitle);
	xlabel('Iteration number');
	ylabel('Error');
	zlabel('Frequency');
end
