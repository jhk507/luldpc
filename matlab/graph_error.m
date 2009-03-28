% $URL$
% $Date$
% $Rev$

function graph_error(filename, titlename, nplot, axis_err)
    % Load the histogram data.
    hist = load(['hist_err_', filename, '.tsv']);
    axis_iters = 0:(size(hist,2)-1);

    % Plot the surface.
    subplot(2,2,nplot);
    surf(axis_iters, axis_err, hist);

    % Make the labels.
    title([titlename, ' error histogram']);
    xlabel('Iteration number');
    ylabel([titlename, ' error']);
    zlabel('Frequency');
end
