% $URL$
% $Date$
% $Rev$

function graph_maxiter(filename, titlename, nplot, axis_err, axis_snr)
    % Load the histogram data.
    hist = load(['hist_maxiter_', filename, '.tsv']);

    % Plot the surface.
    subplot(2,2,nplot);
    surf(axis_err, axis_snr, hist);

    % Make the labels.
    title([titlename, ' maximal iteration error histogram']);
    xlabel([titlename, ' error']);
    ylabel('SNR (dB)');
    zlabel('Frequency');
end
