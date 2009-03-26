% $URL$
% $Date$
% $Rev$

% Load the histogram data.
load histogram.tsv;

% Calculate the histogram dimensions.
niters   = size(histogram,1)-1;
nbuckets = size(histogram,2)-1;

% Extract the histogram matrices.
iters   = histogram(2:(niters+1),1);
buckets = histogram(1,2:(nbuckets+1));
freqs   = histogram(2:(niters+1),2:(nbuckets+1));

% Plot the surface.
surf(buckets, iters, freqs)

% Make the labels.
title('Orthagonal error histogram')
xlabel('Orthagonal error')
ylabel('Iteration number')
zlabel('Frequency')
