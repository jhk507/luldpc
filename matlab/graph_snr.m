% $URL$
% $Date$
% $Rev$

function graph_snr(file, error, nplot)
    % Load the histogram data.
    hist = load(['hist_snr_', file, '.tsv']);

    % Calculate the histogram dimensions.
    niters = size(hist,1)-1;
    nsnrs  = size(hist,2)-1;

    % Extract the histogram matrices.
    iters = hist(2:(niters+1),1);
    snrs  = hist(1,2:(nsnrs+1));
    freqs = hist(2:(niters+1),2:(nsnrs+1));

    % Plot the surface.
    subplot(1,2,nplot)
    surf(snrs, iters, freqs)

    % Make the labels.
    title(['Zero-', error, '-error SNR histogram'])
    xlabel('SNR (dB)')
    ylabel('Iteration number')
    zlabel('Frequency')
end
