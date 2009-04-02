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
		line = fgetl(hdecodes);
		if line == -1
			break;
		end
		axis_decode{d} = line;
		d = d + 1;
	end
	fclose(hdecodes);
	
	global ndecodes;
	ndecodes = length(axis_decode);
	
	% Set all figures to be docked.
	set(0, 'DefaultFigureWindowStyle', 'docked');
	
	% Create the docked figure windows in order.
	for fig = 0:6
		for d = 1:ndecodes
			hfigure = figure(fig*ndecodes + d);
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
		title = ['Block error rate vs. maximum iterations, ', method];
		incfigure(title);
		graph_curves_snr(orthfile, title, axis_snr);

		% Display the SNR/BLER multiple iteration curves
		title = ['Block error rate vs. signal-to-noise ratio, ', method];
		incfigure(title);
		graph_curves_iter(orthfile, title, axis_snr);

		% Display the maxiter/aveiter multiple SNR curves
		title = ['Average iterations to resolution vs. maximum iterations, ', method];
		incfigure(title);
		graph_curves_aveiter(orthfile, title, axis_snr);
		
		% Display the error/iteration/frequency surface
		title = ['Histogram for one SNR, ', method];
		incfigure(title);
		graph_surf_error(orthfile, [title, ', orthagonal error'], axis_err_orth, [-50 2.5 5], 1);
		graph_surf_error(messfile, [title, ', message error'],    axis_err_mess, [-50 0.5 5], 2);
		colorbar;

		% Display the snr/iteration/frequency surface
		title = ['SNR/iteration resolution histogram, ', method];
		incfigure(title);
		graph_surf_snr(orthfile, [title, ', orthagonal error'], axis_snr, 1);
		graph_surf_snr(messfile, [title, ', message error'],    axis_snr, 2);
		colorbar;

		% Display the error/snr/frequency surface
		title = ['Histogram at maximal iteration, ', method];
		incfigure(title);
		graph_surf_maxiter(orthfile, [title, ', orthagonal error'], axis_err_orth, axis_snr, [2    -4 4], 1);
		graph_surf_maxiter(messfile, [title, ', message error'],    axis_err_mess, axis_snr, [0.35 -4 4], 2);
		colorbar;

		% Display the snr/iteration/error/frequency volumetric slices
		title = ['Full volumetric histogram, ', method];
		incfigure(title);
		graph_slice(orthfile, [title, ', orthagonal error'], axis_snr, axis_err_orth, [-300 2.5 2], 1);
		graph_slice(messfile, [title, ', message error'],    axis_snr, axis_err_mess, [-300 0.5 2], 2);
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
