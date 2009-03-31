% $URL$
% $Date$
% $Rev$

function graph_surf_error(filename, titlename, axis_err, cam, nplot)
	% Load the histogram data.
	hist = load(['hist_err_', filename, '.tsv']);
	axis_iters = 0:(size(hist,2)-1);

	% Plot the surface.
	subplot(2,1,nplot);
	surf(axis_iters, axis_err, hist);

	% Make the labels.
	title([titlename, ' error histogram']);
	xlabel('Iteration number');
	ylabel('Error');
	zlabel('Frequency');
	
	axis tight;
	campos(cam);
end
