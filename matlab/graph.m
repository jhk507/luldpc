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
    graph_error('bp_orth',    'BP Orthagonal',    axis_err_orth, [-180 2    4.5], 1);
    graph_error('bp_mess',    'BP Message',       axis_err_mess, [-180 0.35 4.5], 2);
    graph_error('offms_orth', 'OFFMS Orthagonal', axis_err_orth, [-180 2    4.5], 3);
    graph_error('offms_mess', 'OFFMS Message',    axis_err_mess, [-180 0.35 4.5], 4);

    % Display the snr/iteration/frequency surface
    figure(2);
    graph_snr('bp_orth',    'BP Orthagonal',    axis_snr, 1);
    graph_snr('bp_mess',    'BP Message',       axis_snr, 2);
    graph_snr('offms_orth', 'OFFMS Orthagonal', axis_snr, 3);
    graph_snr('offms_mess', 'OFFMS Message',    axis_snr, 4);
    
    % Display the error/snr/frequency surface
    figure(3);
    graph_maxiter('bp_orth',    'BP Orthagonal',    axis_err_orth, axis_snr, [2    -4 4], 1);
    graph_maxiter('bp_mess',    'BP Message',       axis_err_mess, axis_snr, [0.35 -4 4], 2);
    graph_maxiter('offms_orth', 'OFFMS Orthagonal', axis_err_orth, axis_snr, [2    -4 4], 3);
    graph_maxiter('offms_mess', 'OFFMS Message',    axis_err_mess, axis_snr, [0.35 -4 4], 4);
    
    % Display the snr/iteration/error/frequency volumetric slices
    figure(4);
    graph_slice('bp_orth',    'BP Orthagonal',    axis_snr, axis_err_orth, [-450 5 2.5], 1);
    graph_slice('offms_orth', 'OFFMS Orthagonal', axis_snr, axis_err_orth, [-450 5 2.5], 2);
    figure(5);
    graph_slice('bp_mess',    'BP Message',    axis_snr, axis_err_mess, [-900 2.25 4.5], 1);
    graph_slice('offms_mess', 'OFFMS Message', axis_snr, axis_err_mess, [-900 2.25 4.5], 2);
end
