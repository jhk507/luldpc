% $URL$
% $Date$
% $Rev$

function graph_slice(filename, titlename, nplot, axis_snr, axis_err)
    % Load the slice data.
    hist = load(['hist_slice_', filename, '.tsv']);
    
    % Calculate dimensions.
    nsnrs = length(axis_snr);
    niters = size(hist,1)/nsnrs;
    nerrs = length(axis_err);
    
    % Generate the iterator column.
    axis_iter = (0:(niters-1)).';
    
    % Extract the slice matrices.
    vol = zeros(nsnrs, niters, nerrs);
    for s = 1:nsnrs
        vol(s,:,:) = hist( ...
            (1+(s-1)*niters):(s*niters), ...
            1:nerrs);
    end
    
    % Create the stupidly redundant and inexplicably reversed coordinate system.
    [cx, cy, cz] = meshgrid(axis_iter, axis_snr, axis_err);

    % Do the slice plot.
    subplot(2,1,nplot)
    slice(cx, cy, cz, vol, niters-1, axis_snr, 0)

    % Make the labels.
    title([titlename, ' error histogram'])
    xlabel('Iteration number')
    ylabel('SNR (dB)')
    zlabel([titlename, ' error'])
end