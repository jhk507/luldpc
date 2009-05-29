% $URL$
% $Date$
% $Rev$

function graph_curves_iter(hist, axistitle)
	
% Declare the global variables.
	global ndecodes;
	global axis_siter;
	global methods;
	global axis_snr;
	
	% Plot the curves.
	for d = 1:ndecodes
		semilogy(axis_snr, hist(:,d));
		if d < ndecodes
			hold all;
		end
	end
	
	% Make the labels and legend.
	title(axistitle);
	xlabel('SNR (dB)');
	ylabel('BLER');
	legend(methods, 'Location', 'BestOutside');
	
	% Set the axes and grids.
	axis tight;
	grid on;
end
