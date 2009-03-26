% $URL$
% $Date$
% $Rev$

function graph
    clear;
    loadHist('hist_orth.tsv', 'Orthagonal', 1);
    loadHist('hist_mess.tsv', 'Message', 2);
end

function loadHist(file, error, nplot)
    % Load the histogram data.
    hist = load(file);

    % Calculate the histogram dimensions.
    niters   = size(hist,1)-1;
    nbuckets = size(hist,2)-1;

    % Extract the histogram matrices.
    iters   = hist(2:(niters+1),1);
    buckets = hist(1,2:(nbuckets+1));
    freqs   = hist(2:(niters+1),2:(nbuckets+1));

    % Plot the surface.
    subplot(1,2,nplot)
    surf(buckets, iters, freqs)

    % Make the labels.
    title([error, ' error histogram'])
    xlabel([error, ' error'])
    ylabel('Iteration number')
    zlabel('Frequency')
end
