% $URL$
% $Date$
% $Rev$

function graph_scatter(file, error, nplot)
    % Load the scatter data.
    hist = load(['hist_scat_', file, '.tsv']);

    % Extract the scatter matrices.
    snrs    = hist(:,1);
    iters   = hist(:,2);
    buckets = hist(:,3);
    sizes   = hist(:,4)*40 + 1;
    colours = hist(:,4);

    % Do the scatterplot.
    subplot(2,1,nplot)
    scatter3(snrs, iters, buckets, sizes, colours, 'filled')

    % Make the labels.
    title([error, ' error histogram'])
    xlabel('SNR (dB)')
    ylabel('Iteration number')
    zlabel([error, ' error'])
end
