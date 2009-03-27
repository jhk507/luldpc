% $URL$
% $Date$
% $Rev$

function graph_slice(filename, titlename, nplot, axis_snr, axis_err)
    % Load the slice data.
    hist = load(['hist_slice_', filename, '.tsv']);

    % Extract the slice matrices.

    % Do the slice plot.
    subplot(2,1,nplot)

    % Make the labels.
    title([titlename, ' error histogram'])
    xlabel('SNR (dB)')
    ylabel('Iteration number')
    zlabel([titlename, ' error'])
end
