% $URL$
% $Date$
% $Rev$

function graph_curves_aveiter_snr(aves, axistitle)
	% Declare the global variables.
	global ndecodes axis_iter methods;
	
	% Plot the curves.
	for d = 1:ndecodes
		plot(axis_iter(1:10:end), aves(d,1:10:end));
		if d < ndecodes
			hold all;
		end
	end
	
	% Make the labels and legend.
	title(axistitle);
	xlabel('Maximum Iteration Number');
	ylabel('Average Iteration Number');;
	set(legend(methods, 'Location', 'BestOutside'), 'Interpreter','none');
	
	% Set the axes and grids.
	axis tight;
	grid on;
end



