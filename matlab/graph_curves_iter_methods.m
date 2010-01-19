% $URL$
% $Date$
% $Rev$

function graph_curves_iter_methods(hist, axistitle)
	% Declare the global variables.
	global nsiters nsnrs axis_siter axis_siter_text axis_snr;
	
	% Plot the curves.
	for i = 1:nsiters
		semilogy(axis_snr, hist(:,axis_siter(i)));
		if i < nsiters
			hold all;
		end
	end
	
	% Make the labels and legend.
	title(axistitle);
	xlabel('Eb/N0 (dB)');
	ylabel('BLER');
	set(legend(axis_siter_text, 'Location', 'BestOutside'), 'Interpreter','none');
	
	% Set the axes and grids.
	axis tight;
	grid on;
end
