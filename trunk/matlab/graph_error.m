% $URL$
% $Date$
% $Rev$

function graph_error(filename, titlename, nplot, axis_err)
    % Load the histogram data.
    hist = load(['hist_err_', filename, '.tsv']);
    axis_iters = 1:size(hist,1);

    % Plot the surface.
    subplot(2,2,nplot)
    surf(axis_err, axis_iters, hist)

    % Make the labels.
    title([titlename, ' error histogram'])
    xlabel([titlename, ' error'])
    ylabel('Iteration number')
    zlabel('Frequency')
end
