% $URL$
% $Date$
% $Rev$

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the histogram data.
load hist_orth.tsv;

% Calculate the histogram dimensions.
niters   = size(hist_orth,1)-1;
nbuckets = size(hist_orth,2)-1;

% Extract the histogram matrices.
iters   = hist_orth(2:(niters+1),1);
buckets = hist_orth(1,2:(nbuckets+1));
freqs   = hist_orth(2:(niters+1),2:(nbuckets+1));

% Plot the surface.
subplot(1,2,1)
surf(buckets, iters, freqs)

% Make the labels.
title('Orthagonal error histogram')
xlabel('Orthagonal error')
ylabel('Iteration number')
zlabel('Frequency')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the histogram data.
load hist_mess.tsv;

% Calculate the histogram dimensions.
niters   = size(hist_mess,1)-1;
nbuckets = size(hist_mess,2)-1;

% Extract the histogram matrices.
iters   = hist_mess(2:(niters+1),1);
buckets = hist_mess(1,2:(nbuckets+1));
freqs   = hist_mess(2:(niters+1),2:(nbuckets+1));

% Plot the surface.
subplot(1,2,2)
surf(buckets, iters, freqs)

% Make the labels.
title('Message error histogram')
xlabel('Message error')
ylabel('Iteration number')
zlabel('Frequency')
