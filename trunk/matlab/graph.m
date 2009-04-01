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
		line = {fgetl(hdecodes)}
		if line{1} == -1
			break;
		end
		axis_decode(d) = line;
		d = d + 1;
	end
	fclose(hdecodes);
	
	global ndecodes;
	ndecodes = length(axis_decode);
	
	% Loop through the decode methods.
	for d = 1:ndecodes
		orthfile = [axis_decode{d}, '_orth'];
		messfile = [axis_decode{d}, '_mess'];
		title = upper(axis_decode{d});
		
		global f;
		f = d;
		
		% Display the iteration/BLER multiple SNR curves
		incfigure;
		graph_curves_snr(orthfile, title, axis_snr);

		% Display the SNR/BLER multiple iteration curves
		incfigure;
		graph_curves_iter(orthfile, title, axis_snr);

		% Display the maxiter/aveiter multiple SNR curves
		incfigure;
		graph_curves_aveiter(orthfile, title, axis_snr);
		
		orthtitle = [title, ' Orthagonal'];
		messtitle = [title, ' Message'];
		
		% Display the error/iteration/frequency surface
		incfigure;
		graph_surf_error(orthfile, orthtitle, axis_err_orth, [-180 2    4.5], 1);
		graph_surf_error(messfile, messtitle, axis_err_mess, [-180 0.35 4.5], 2);

		% Display the snr/iteration/frequency surface
		incfigure;
		graph_surf_snr(orthfile, orthtitle, axis_snr, 1);
		graph_surf_snr(messfile, messtitle, axis_snr, 2);

		% Display the error/snr/frequency surface
		incfigure;
		graph_surf_maxiter(orthfile, orthtitle, axis_err_orth, axis_snr, [2	   -4 4], 1);
		graph_surf_maxiter(messfile, messtitle, axis_err_mess, axis_snr, [0.35 -4 4], 2);

		% Display the snr/iteration/error/frequency volumetric slices
		incfigure;
		graph_slice(orthfile, orthtitle, axis_snr, axis_err_orth, [-450 5    2.5], 1);
		graph_slice(messfile, messtitle, axis_snr, axis_err_mess, [-900 2.25 4.5], 2);
	end
end

% Increment the figure number and set the window to docked mode.
function incfigure% (methodtitle, 
	global f;
	hfigure = figure(f);
	clf;
	set(hfigure, 'WindowStyle', 'docked');
	
	axes = get(hfigure, 'CurrentAxes');
%	title = get(axes, 'Title');
%	set(hfigure, 'Name', title);

	global ndecodes;
	f = f+ndecodes;
end
