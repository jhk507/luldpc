% $URL$
% $Date$
% $Rev$

% Do all graphs.
function graph
	% Clear previous workspace data
	clear;
	
	% Load the axes.
	axis_snr = load('axis_snr.tsv');
	axis_err_orth = load('axis_err_orth.tsv');
	axis_err_mess = load('axis_err_mess.tsv');
	
	% Load the decode method names.
	hdecodes = fopen('axis_decode.tsv', 'r');
	d = 1;
	while 1
		line = {fgetl(hdecodes)};
		if line{1} == -1
			break;
		end
		axis_decode(d) = line;
		d = d + 1;
	end
	fclose(hdecodes);
	
	global ndecodes;
	ndecodes = length(axis_decode);
	
	% Create the docked figure windows in order.
	for fig = 0:6
		for d = 1:ndecodes
			hfigure = figure(fig*ndecodes + d);
			set(hfigure, 'WindowStyle', 'docked');
			clf;
		end
	end

	global f;
	
	% Loop through the decode methods.
	for d = 1:ndecodes
		orthfile = [axis_decode{d}, '_orth'];
		messfile = [axis_decode{d}, '_mess'];
		method = upper(axis_decode{d});
		
		f = d;
		
		% Display the iteration/BLER multiple SNR curves
		title = [method, ' block error rate vs. maximum iterations'];
		incfigure(title);
		graph_curves_snr(orthfile, title, axis_snr);

		% Display the SNR/BLER multiple iteration curves
		title = [method, ' block error rate vs. signal-to-noise ratio'];
		incfigure(title);
		graph_curves_iter(orthfile, title, axis_snr);

		% Display the maxiter/aveiter multiple SNR curves
		title = [method, ' average iterations to resolution vs. maximum iterations'];
		incfigure(title);
		graph_curves_aveiter(orthfile, title, axis_snr);
		
		% Display the error/iteration/frequency surface
		title = [method, ' error histogram'];
		incfigure(title);
		graph_surf_error(orthfile, ['Orthagonal ', title], axis_err_orth, [-50 2.5 5], 1);
		graph_surf_error(messfile, ['Message ',    title], axis_err_mess, [-50 0.5 5], 2);
		colorbar;

		% Display the snr/iteration/frequency surface
		title = [method, ' zero-error SNR histogram'];
		incfigure(title);
		graph_surf_snr(orthfile, ['Orthagonal ', title], axis_snr, 1);
		graph_surf_snr(messfile, ['Message ',    title], axis_snr, 2);
		colorbar;

		% Display the error/snr/frequency surface
		title = [method, ' maximal iteration error histogram'];
		incfigure(title);
		graph_surf_maxiter(orthfile, ['Orthagonal ', title], axis_err_orth, axis_snr, [2    -4 4], 1);
		graph_surf_maxiter(messfile, ['Message ',    title], axis_err_mess, axis_snr, [0.35 -4 4], 2);
		colorbar;

		% Display the snr/iteration/error/frequency volumetric slices
		title = [method, ' error histogram'];
		incfigure(title);
		graph_slice(orthfile, ['Orthagonal ', title], axis_snr, axis_err_orth, [-300 2.5 2], 1);
		graph_slice(messfile, ['Message ',    title], axis_snr, axis_err_mess, [-300 0.5 2], 2);
		colorbar;
	end
end

% Increment the figure number and set the window to docked mode.
function incfigure(figtitle)
	global f;
	hfigure = figure(f);
	set(hfigure, 'Name', figtitle);

	global ndecodes;
	f = f+ndecodes;
end
