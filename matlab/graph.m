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
	global nsiters;	nsiters = 2;
	global niters; 	niters  = length(axis_iter);
	global nsnrs;   nsnrs   = length(axis_snr);

	% Create the iteration legend text.
	global axis_siter;      axis_siter = zeros(nsiters, 1);
	global axis_siter_text; axis_siter_text = cell(nsiters, 1);
	for i = 1:nsiters
		axis_siter(i) = floor(niters/10*(9*(i-1)/(nsiters-1)+1));
		axis_siter_text{i} = ['i=', num2str(axis_siter(i))];
	end
	
	% Create the SNR text for the axis.
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
		axis_decode{d} = line;
		d = d + 1;
	end
	fclose(hdecodes);
	
	% Load the decode-related information and construct file names.
	global ndecodes;	ndecodes = length(axis_decode);

	global orthfiles;	orthfiles = cell(ndecodes, 1);
	global messfiles;	messfiles = cell(ndecodes, 1);
	global methods;		methods = cell(ndecodes, 1);
	
	for d = 1:ndecodes
		orthfiles{d} = [axis_decode{d}, '_orth'];
		messfiles{d} = [axis_decode{d}, '_mess'];
		methods{d} = upper(axis_decode{d});
	end
	
	% Calculate the 'average iterations' data
	hist_aveiter = zeros(nsnrs, ndecodes, niters);
	hist_aveiter_alt = zeros(nsnrs, ndecodes, niters);
	hist_maxiter = zeros(nsnrs, ndecodes, niters);
	hist_bleriter = zeros(nsnrs, ndecodes, niters);
	
	for d = 1:ndecodes
		hist = load(['hist_snr_', orthfiles{d}, '.tsv']);
		
		% Convert from a cumulative histogram to a moving weighted average
		wsum = zeros(nsnrs, 1);
		for i = 1:niters
			weight = hist(:,i);
			if i > 1
				weight = weight - hist(:,i-1);
			end
			weight = weight*i;
			wsum = wsum + weight;
			weight = wsum;
			for s = 1:nsnrs
				if hist(s,i) ~= 0
					weight(s) = weight(s)./hist(s,i);
				end
				
			end
			hist_aveiter(:,d,i) = weight;
		end
		hist_maxiter(:,d,:) = 1 - load(['hist_snr_', orthfiles{d}, '.tsv']);
		hist_bleriter(:,d,:) = 1 - load(['hist_snr_', orthfiles{d}, '.tsv']);
		
		hist_aveiter_alt(:,d,:) = load(['hist_iter_', axis_decode{d}, '.tsv']);
	end
	
	% Set all figures to be docked.
	set(0, 'DefaultFigureWindowStyle', 'docked');
	
	nfigures = 4*ndecodes;
	
	global f;
	
	% Create the docked figure windows in order.
	for f = 1:nfigures
		figure(f);
	end
	
	% Loop through the decode methods.
	for d = 1:ndecodes
		f = d;
		
		% Display the iteration/BLER multiple SNR curves
		title = ['Block Error Rate vs. Maximum Iterations, ', methods{d}];
		incfigure(title, ndecodes);
		graph_curves_snr(squeeze(hist_maxiter(:,d,:)), title);

		% Display the SNR/BLER multiple iteration curves
		title = ['Block Error Rate vs. Signal-to-noise Ratio, ', methods{d}];
		incfigure(title, ndecodes);
		graph_curves_iter_methods(squeeze(hist_bleriter(:,d,:)), title);

		% Display the maxiter/aveiter multiple SNR curves
		title = ['Average Iteration Number vs. Maximum Iterations Number, ', methods{d}];
		incfigure(title, ndecodes);
		graph_curves_aveiter(squeeze(hist_aveiter(:,d,:)), title);
		
		% Display the error/iteration/frequency surface
		%title = ['Histogram for one SNR, ', methods{d}];
		%incfigure(title, ndecodes);
		%graph_surf_error(orthfiles{d}, [title, ', orthagonal error'], axis_err_orth, 1);
		%graph_surf_error(messfiles{d}, [title, ', message error'],    axis_err_mess, 2);
		%colorbar;

		% Display the snr/iteration/frequency surface
		%title = ['SNR/iteration resolution histogram, ', methods{d}];
		%incfigure(title, ndecodes);
		%graph_surf_snr(orthfiles{d}, [title, ', orthagonal error'], 1);
		%graph_surf_snr(messfiles{d}, [title, ', message error'],    2);
		%colorbar;

		% Display the error/snr/frequency surface
		%title = ['Histogram at maximal iteration, ', methods{d}];
		%incfigure(title, ndecodes);
		%graph_surf_maxiter(orthfiles{d}, [title, ', orthagonal error'], axis_err_orth, [2    -4 4], 1);
		%graph_surf_maxiter(messfiles{d}, [title, ', message error'],    axis_err_mess, [0.35 -4 4], 2);
		%colorbar;
		
		% Display the snr/iteration/error averaged histogram surface
		%title = ['Weighted mean histogram, ', methods{d}];
		%incfigure(title, ndecodes);
		%graph_surf_histave(orthfiles{d}, [title, ', orthagonal error'], axis_err_orth, 1);
		%graph_surf_histave(messfiles{d}, [title, ', message error'],    axis_err_mess, 2);

		%Display the snr/iteration/error/frequency volumetric slices
		title = ['Full volumetric histogram, ', methods{d}];
		incfigure(title, ndecodes);
		graph_slice(orthfiles{d}, [title, ', orthagonal error'], axis_err_orth, [-300 2.5 2], 1);
		graph_slice(messfiles{d}, [title, ', message error'],    axis_err_mess, [-300 0.5 2], 2);
		colorbar;
	end

	f = nfigures+1;
	
	displaysnrs = [ 1.5 1.8 ];
	
	for dsnri = 1:length(displaysnrs)
		snr = displaysnrs(dsnri);
		for snri = 1:nsnrs
			if axis_snr(snri) == snr
				% Display the iteration/BLER multiple SNR curves
				title = ['Block Error Rate vs. Maximum Iterations, ', axis_snr_text{snri}];
				incfigure(title, 1);
				graph_curves_snr_methods(reshape(hist_maxiter(snri,:,:),ndecodes,niters), title);	
				
				% Display the maxiter/aveiter multiple SNR curves
				title = ['Average Iteration Number vs. Maximum Iteration Number, ', axis_snr_text{snri}];
				incfigure(title, 1);
				graph_curves_aveiter_snr(reshape(hist_aveiter(snri,:,:),ndecodes,niters), title);
				
				% Display the maxiter/aveitr multiple SNR curves
				% (alternate method)
				title = ['Average Iteration Number vs. Maximum Iteration Number (alternate method), ', axis_snr_text{snri}];
				incfigure(title, 1);
				graph_curves_aveiter_snr(reshape(hist_aveiter_alt(snri,:,:),ndecodes,niters), title);
				
				break;
			end
		end
	end

	for i= 1:nsiters
		iter= axis_siter(i);
		
		title = ['Block Error Rate vs. Signal-to-noise Ratio, i=', num2str(iter)];
		incfigure(title, 1);
		graph_curves_iter(squeeze(hist_bleriter(:,:,iter)), title);	
	end
	
	% Display the performance histogram
	title = 'Single-iteration decode performance histogram';
	incfigure(title, 1);
	graph_curves_perf(title);
end

% Increment the figure number and set the window to docked mode.
function incfigure(figtitle, increment)
	global f;
	hfigure = figure(f);
	clf;
	set(hfigure, 'Name', figtitle);

	f = f+increment;
end
