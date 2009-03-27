% $URL$
% $Date$
% $Rev$

function graph_slice(filename, titlename, nplot, axis_snr, axis_err)
    % Load the slice data.
    hist = load(['hist_slice_', filename, '.tsv']);
    
    % Calculate dimensions.
    nsnrs = size(axis_snr,1);
    niters = size(hist,1)/nsnrs;
    nerrs = size(axis_err,1);
    
    % Generate the iterator column.
    axis_iters = (0:(niters-1)).';
    
    % Extract the slice matrices.
    vol = zeros(nsnrs, niters, nerrs);
    for s = 1:nsnrs
        vol(s,:,:) = hist( ...
            (1+(s-1)*niters):(s*niters), ...
            1:nerrs);
    end

    % Do the slice plot.
    subplot(2,1,nplot)
    slice(vol, axis_snr, niters-1, 0);

    % Make the labels.
    title([titlename, ' error histogram'])
    xlabel('SNR (dB)')
    ylabel('Iteration number')
    zlabel([titlename, ' error'])
end
