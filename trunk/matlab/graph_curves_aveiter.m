% $URL$
% $Date$
% $Rev$

function graph_curves_aveiter(aves, axistitle)
	% Declare the global variables.
	global nsnrs;
	global axis_iter;
	global axis_snr_text;
	
	% Plot the curves.
	for s = 1:nsnrs
		plot(axis_iter, aves(s,:));
		if s < nsnrs
			hold all;
		end
	end
	
	% Make the labels and legend.
	title(axistitle);
	xlabel('Maximum iteration number');
	ylabel('Average iteration number');
	set(legend(axis_snr_text, 'Location', 'BestOutside'), 'Interpreter','none');

	% Set the axes and grids.
	axis tight;
	grid on;
end
