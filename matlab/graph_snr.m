% $URL$
% $Date$
% $Rev$

function graph_snr(filename, titlename, nplot, axis_snr)
    % Load the histogram data.
    hist = load(['hist_snr_', filename, '.tsv']);
    axis_iters = 0:(size(hist,2)-1);

    % Plot the surface.
    subplot(2,2,nplot);
    surf(axis_snr, axis_iters, hist);

    % Make the labels.
    title([titlename, ' zero-error SNR histogram']);
    xlabel('Iteration number');
    ylabel('SNR (dB)');
    zlabel('Frequency');
end
