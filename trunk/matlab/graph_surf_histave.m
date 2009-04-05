% $URL$
% $Date$
% $Rev$

function graph_surf_histave(filename, axistitle, axis_err, nplot)
	% Load the slice data.
	hist = load(['hist_slice_', filename, '.tsv']);
	
	% Declare the global variables.
	global niters;
	global nsnrs;
	global axis_iter;
	global axis_snr;
	
	nerrs = length(axis_err);
	
	% Compute the weighted average.
	avehist = zeros(niters, nsnrs);
	sslice = zeros(nerrs, niters);
	for s = 1:nsnrs
		sslice = hist((1+(s-1)*nerrs):(s*nerrs), 1:niters);
		avehist(:,s) = (axis_err.' * sslice) ./ sum(sslice, 1);
	end
	
	% Plot the surface.
	subplot(1,2,nplot);
	surf(axis_snr, axis_iter, avehist, 'MeshStyle', 'column');
	
	% Set up the axis properties.
	set(gca, 'XDir', 'reverse', 'YDir', 'reverse', 'YScale', 'log');
	axis tight;

	% Make the labels and colorbar.
	title(axistitle);
	xlabel('SNR (dB)');
	ylabel('Iteration number');
	zlabel('Error');
	colorbar;
end
