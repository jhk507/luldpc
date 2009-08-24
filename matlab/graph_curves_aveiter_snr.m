% $URL$
% $Date$
% $Rev$

function graph_curves_aveiter_snr(aves, axistitle)
	% Declare the global variables.
	global ndecodes;
	global axis_iter;
	global methods;
	
	% Plot the curves.
	for d = 1:ndecodes
		plot(axis_iter, aves(d,:));
		if d < ndecodes
			hold all;
		end
	end
	
	% Make the labels and legend.
	title(axistitle);
	xlabel('Maximum iterations');
	ylabel('Average iterations to resolution');
	set(legend(methods, 'Location', 'BestOutside'), 'Interpreter','none');
	
	% Set the axes and grids.
	axis tight;
	grid on;
end