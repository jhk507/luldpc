% $URL$
% $Date$
% $Rev$

function graph_slice(filename, axistitle, axis_err, cam, nplot)
	% Load the slice data.
	hist = load(['hist_slice_', filename, '.tsv']);
	
	% Declare the global variables.
	global niters;
	global nsnrs;
	global axis_iter;
	global axis_snr;
	
	nerrs = length(axis_err);
	
	% Extract the slice matrices.
	vol = zeros(nerrs, niters, nsnrs);
	for s = 1:nsnrs
		vol(:,:,s) = hist((1+(s-1)*nerrs):(s*nerrs), 1:niters);
	end
	
	% Create the redundant and inexplicably reversed coordinate system.
	[cx, cy, cz] = meshgrid(axis_iter, axis_err, axis_snr);

	% Do the slice plot.
	subplot(1,2,nplot);
	hslice = slice(cx, cy, cz, vol, niters-1, 0, axis_snr);
	
	% Set up the axis properties.
	set(gca, 'XScale', 'log', 'YDir', 'reverse');
	axis tight;

	% Make the labels.
	title(axistitle);
	xlabel('Iteration number');
	ylabel('Error');
	zlabel('SNR (dB)');

	% Remove some lines on the slices.
	for hi = 1:length(hslice)
		dims = size(get(hslice(hi), 'XData'));
		if dims == [niters nsnrs]
			set(hslice(hi), 'MeshStyle', 'column');
		elseif dims == [nerrs nsnrs]
			set(hslice(hi), 'MeshStyle', 'column');
		elseif dims == [nerrs niters]
			set(hslice(hi), 'LineStyle', 'none');
		end
	end

	% Set up the camera.
	campos(cam);
end
