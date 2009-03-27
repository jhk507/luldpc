% $URL$
% $Date$
% $Rev$

function graph
    % Clear previous workspace data
    clear;
    
    % Load the axes.
    axis_snr = load('axis_snr.tsv');
    axis_err_orth = load('axis_err_orth.tsv');
    axis_err_mess = load('axis_err_mess.tsv');

    % Display the error/iteration/frequency surface
    figure(1);
    graph_error('bp_orth',    'BP Orthagonal',    1, axis_err_orth);
    graph_error('bp_mess',    'BP Message',       2, axis_err_mess);
    graph_error('offms_orth', 'OFFMS Orthagonal', 3, axis_err_orth);
    graph_error('offms_mess', 'OFFMS Message',    4, axis_err_mess);

    % Display the snr/iteration/frequency surface
    figure(2);
    graph_snr('bp_orth',    'BP Orthagonal',    1, axis_snr);
    graph_snr('bp_mess',    'BP Message',       2, axis_snr);
    graph_snr('offms_orth', 'OFFMS Orthagonal', 3, axis_snr);
    graph_snr('offms_mess', 'OFFMS Message',    4, axis_snr);
    
    % Display the snr/iteration/error/frequency surface
    figure(3);
    graph_scatter('bp_orth',    'BP Orthagonal',    1, axis_snr, axis_err_orth);
    graph_scatter('offms_orth', 'OFFMS Orthagonal', 2, axis_snr, axis_err_orth);
    figure(4)
    graph_scatter('bp_mess',    'BP Message',    1, axis_snr, axis_err_mess);
    graph_scatter('offms_mess', 'OFFMS Message', 2, axis_snr, axis_err_mess);
end
