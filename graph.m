% $URL$
% $Date$
% $Rev$

% Load the histogram data.
load histogram.tsv;

% Calculate the histogram dimensions.
niters   = size(hist,1)-1;
nbuckets = size(hist,2)-1;

% Extract the histogram matrices.
iters   = hist(2:(niters+1),1);
buckets = hist(1,2:(nbuckets+1));
freqs   = hist(2:(niters+1),2:(nbuckets+1));

% Plot the surface.
surf(buckets, iters, freqs)
