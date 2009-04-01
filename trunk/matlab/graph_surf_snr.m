% $URL$
% $Date$
% $Rev$

function graph_surf_snr(filename, titlename, axis_snr, nplot)
	% Load the histogram data.
	hist = load(['hist_snr_', filename, '.tsv']);
	axis_iters = 0:(size(hist,2)-1);

	% Plot the surface.
	subplot(2,1,nplot);
	hsurf = surf(axis_iters, axis_snr, hist);

	% Make the labels.
	title([titlename, ' zero-error SNR histogram']);
	xlabel('Iteration number');
	ylabel('SNR (dB)');
	zlabel('Frequency');

	axis tight;
	set(hsurf, 'MeshStyle', 'row');
end
