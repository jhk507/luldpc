% $URL$
% $Date$
% $Rev$

function graph_curves_snr_methods(hist, axistitle)
	% Declare the global variables.
	global ndecodes axis_iter methods;
	
	% Plot the curves.
	for d = 1:ndecodes
		semilogy(axis_iter(1:10:end), hist(d,1:10:end));
		if d < ndecodes
			hold all;
		end
	end
	
	% Make the labels and legend.
	title(axistitle);
	xlabel('Iteration number');
	ylabel('BLER');
	set(legend(methods, 'Location', 'BestOutside'), 'Interpreter','none');
	
	% Set the axes and grids.
	axis tight;
	grid on;
end