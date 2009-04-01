% $URL$
% $Date$
% $Rev$

function graph_surf_maxiter(filename, titlename, axis_err, axis_snr, cam, nplot)
	% Load the histogram data.
	hist = load(['hist_maxiter_', filename, '.tsv']);

	% Plot the surface.
	subplot(2,1,nplot);
	hsurf = surf(axis_err, axis_snr, hist);

	% Make the labels.
	title([titlename, ' maximal iteration error histogram']);
	xlabel('Error');
	ylabel('SNR (dB)');
	zlabel('Frequency');
	
	axis tight;
	campos(cam);
	set(hsurf, 'MeshStyle', 'both');
end
