% $URL$
% $Date$
% $Rev$

function graph
    % Clear previous workspace data
    clear;

    % Display the error/iteration/frequency surface
    figure(1);
    graph_error('orth', 'Orthagonal', 1);
    graph_error('mess', 'Message', 2);

    % Display the snr/iteration/frequency surface
    figure(2);
    graph_snr('orth', 'Orthagonal', 1);
    graph_snr('mess', 'Message', 2);
    
    % Display the snr/iteration/error/frequency surface
    figure(3);
    graph_scatter('orth', 'Orthagonal', 1);
    graph_scatter('mess', 'Message', 2);
end
