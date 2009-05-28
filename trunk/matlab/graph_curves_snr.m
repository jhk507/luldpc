% $URL$
% $Date$
% $Rev$

function graph_curves_snr(hist, axistitle)
	% Declare the global variables.
	global nsnrs;
	global axis_iter;
	global axis_snr_text;
	
	% Plot the curves.
	for s = 1:nsnrs
		semilogy(axis_iter, hist(s,:));
		if s < nsnrs
			hold all;
		end
	end
	
	% Make the labels and legend.
	title(axistitle);
	xlabel('Iteration number');
	ylabel('BLER');
	legend(axis_snr_text, 'Location', 'BestOutside');
	
	% Set the axes and grids.
	axis tight;
	grid on;
end
