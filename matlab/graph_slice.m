% $URL$
% $Date$
% $Rev$

function graph_slice(filename, titlename, axis_snr, axis_err, cam, nplot)
	% Load the slice data.
	hist = load(['hist_slice_', filename, '.tsv']);
	
	% Calculate dimensions.
	nerrs = length(axis_err);
	nsnrs = length(axis_snr);
	niters = size(hist,2);
	
	% Generate the iterator column.
	axis_iter = (0:(niters-1)).';
	
	% Extract the slice matrices.
	vol = zeros(nerrs, niters, nsnrs);
	for s = 1:nsnrs
		vol(:,:,s) = hist( ...
			(1+(s-1)*nerrs):(s*nerrs), ...
			1:niters);
	end
	
	% Create the stupidly redundant and inexplicably reversed coordinate system.
	[cx, cy, cz] = meshgrid(axis_iter, axis_err, axis_snr);

	% Do the slice plot.
	subplot(1,2,nplot);
	hslice = slice(cx, cy, cz, vol, niters-1, 0, axis_snr);

	% Make the labels.
	title([titlename, ' error histogram']);
	xlabel('Iteration number');
	ylabel('Error');
	zlabel('SNR (dB)');

	% Set the axis style.
	axis tight;
	% Set the camera position.
	campos(cam);
	% Remove lines on the SNR slices.
	for hi = 1:length(hslice)
		if size(get(hslice(hi), 'XData')) == [nerrs niters]
			set(hslice(hi), 'LineStyle', 'none');
		end
	end
end
