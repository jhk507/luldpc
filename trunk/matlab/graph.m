% $URL$
% $Date$
% $Rev$

function graph
    % Clear previous workspace data
    clear;

    % Display the error/iteration/frequency surface
    figure(1);
    graph_error('bp_orth',    'BP Orthagonal',    1);
    graph_error('bp_mess',    'BP Message',       2);
    graph_error('offms_orth', 'OFFMS Orthagonal', 3);
    graph_error('offms_mess', 'OFFMS Message',    4);

    % Display the snr/iteration/frequency surface
    figure(2);
    graph_snr('bp_orth',    'BP Orthagonal',    1);
    graph_snr('bp_mess',    'BP Message',       2);
    graph_snr('offms_orth', 'OFFMS Orthagonal', 3);
    graph_snr('offms_mess', 'OFFMS Message',    4);
    
    % Display the snr/iteration/error/frequency surface
    figure(3);
    graph_scatter('bp_orth',    'BP Orthagonal',    1);
    graph_scatter('offms_orth', 'OFFMS Orthagonal', 2);
    figure(4)
    graph_scatter('bp_mess',    'BP Message',    1);
    graph_scatter('offms_mess', 'OFFMS Message', 2);
end
