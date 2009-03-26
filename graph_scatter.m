% $URL$
% $Date$
% $Rev$

function graph_scatter
    clear;
    loadScatter('hist_scat_orth.tsv', 'Orthagonal', 1);
    loadScatter('hist_scat_mess.tsv', 'Message', 2);
end

function loadScatter(file, error, nplot)
    % Load the scatter data.
    scat = load(file);

    % Extract the scatter matrices.
    snrs    = scat(:,1);
    iters   = scat(:,2);
    buckets = scat(:,3);
    sizes   = scat(:,4)*80 + 1;
    colours = scat(:,5);

    % Do the scatterplot.
    subplot(2,1,nplot)
    scatter3(snrs, iters, buckets, sizes, colours, 'filled')

    % Make the labels.
    title([error, ' error histogram'])
    xlabel('SNR (dB)')
    ylabel('Iteration number')
    zlabel([error, ' error'])
end
