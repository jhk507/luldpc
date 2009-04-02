% $URL$
% $Date$
% $Rev$

% Do all graphs.
function graph
	% Clear previous workspace data
	clear;
	
	% Load the axes.
	global axis_iter;
	global axis_snr;
	axis_iter     = load('axis_iter.tsv');
	axis_snr      = load('axis_snr.tsv');
	axis_err_orth = load('axis_err_orth.tsv');
	axis_err_mess = load('axis_err_mess.tsv');
	
	% Calculate the axis lengths.
	global nsiters;	nsiters = 10;
	global niters; 	niters  = length(axis_iter);
	global nsnrs;   nsnrs   = length(axis_snr);

	% Create the legend text.
	global axis_siter;      axis_siter = zeros(nsiters, 1);
	global axis_siter_text; axis_siter_text = cell(nsiters, 1);
	for i = 1:nsiters
		axis_siter(i) = max(1, floor(i/nsiters*niters));
		axis_siter_text{i} = ['i=', num2str(axis_siter(i))];
	end
	
	global axis_snr_text; axis_snr_text = cell(nsnrs, 1);
	for s = 1:nsnrs
		axis_snr_text{s} = [num2str(axis_snr(s)),'dB'];
	end

	% Load the decode method names.
	hdecodes = fopen('axis_decode.tsv', 'r');
	global axis_decode;
	d = 1;
	while 1
		line = fgetl(hdecodes);
		if line == -1
			break;
		end
		axis_decode{d} = upper(line);
		d = d + 1;
	end
	fclose(hdecodes);
	
	global ndecodes;
	ndecodes = length(axis_decode);
	
	% Set all figures to be docked.
	set(0, 'DefaultFigureWindowStyle', 'docked');
	
	nfigures = 7*ndecodes;
	
	global f;
	
	% Create the docked figure windows in order.
	for f = 1:nfigures
		figure(f);
		clf;
	end

	% Loop through the decode methods.
	for d = 1:ndecodes
		orthfile = [axis_decode{d}, '_orth'];
		messfile = [axis_decode{d}, '_mess'];
		method = axis_decode{d};
		
		f = d;
		
		% Display the iteration/BLER multiple SNR curves
		title = ['Block error rate vs. maximum iterations, ', method];
		incfigure(title);
		graph_curves_snr(orthfile, title);

		% Display the SNR/BLER multiple iteration curves
		title = ['Block error rate vs. signal-to-noise ratio, ', method];
		incfigure(title);
		graph_curves_iter(orthfile, title);

		% Display the maxiter/aveiter multiple SNR curves
		title = ['Average iterations to resolution vs. maximum iterations, ', method];
		incfigure(title);
		graph_curves_aveiter(orthfile, title);
		
		% Display the error/iteration/frequency surface
		title = ['Histogram for one SNR, ', method];
		incfigure(title);
		graph_surf_error(orthfile, [title, ', orthagonal error'], axis_err_orth, 1);
		graph_surf_error(messfile, [title, ', message error'],    axis_err_mess, 2);
		colorbar;

		% Display the snr/iteration/frequency surface
		title = ['SNR/iteration resolution histogram, ', method];
		incfigure(title);
		graph_surf_snr(orthfile, [title, ', orthagonal error'], 1);
		graph_surf_snr(messfile, [title, ', message error'],    2);
		colorbar;

		% Display the error/snr/frequency surface
		title = ['Histogram at maximal iteration, ', method];
		incfigure(title);
		graph_surf_maxiter(orthfile, [title, ', orthagonal error'], axis_err_orth, [2    -4 4], 1);
		graph_surf_maxiter(messfile, [title, ', message error'],    axis_err_mess, [0.35 -4 4], 2);
		colorbar;

		% Display the snr/iteration/error/frequency volumetric slices
		title = ['Full volumetric histogram, ', method];
		incfigure(title);
		graph_slice(orthfile, [title, ', orthagonal error'], axis_err_orth, [-300 2.5 2], 1);
		graph_slice(messfile, [title, ', message error'],    axis_err_mess, [-300 0.5 2], 2);
		colorbar;
	end
	
	% Display the performance histogram
	title = 'Performance histogram';
	f = nfigures+1;
	incfigure(title);
	graph_curves_perf(title);	
end

% Increment the figure number and set the window to docked mode.
function incfigure(figtitle)
	global f;
	hfigure = figure(f);
	set(hfigure, 'Name', figtitle);

	global ndecodes;
	f = f+ndecodes;
end
