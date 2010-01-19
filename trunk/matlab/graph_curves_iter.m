% $URL$
% $Date$
% $Rev$

function graph_curves_iter(hist, axistitle)
	% Declare the global variables.
	global ndecodes axis_siter  methods axis_snr;
	
	% Plot the curves.
	for d = 1:ndecodes
		semilogy(axis_snr, hist(:,d));
		if d < ndecodes
			hold all;
		end
	end
	
	% Make the labels and legend.
	title(axistitle);
	xlabel('Eb/N0 (dB)');
	ylabel('BLER');
	set(legend(methods, 'Location', 'BestOutside'), 'Interpreter','none');
	
	% Set the axes and grids.
	axis tight;
	grid on;
end
