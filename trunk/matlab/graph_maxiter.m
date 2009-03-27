% $URL$
% $Date$
% $Rev$

function graph_maxiter(filename, titlename, nplot, axis_err, axis_snr)
    % Load the histogram data.
    hist = load(['hist_maxiter_', filename, '.tsv']);

    % Plot the surface.
    subplot(2,2,nplot)
    surf(axis_snr, axis_err, hist)

    % Make the labels.
    title([titlename, ' maximal iteration error histogram'])
    xlabel('SNR (dB)')
    ylabel([titlename, ' error'])
    zlabel('Frequency')
end
