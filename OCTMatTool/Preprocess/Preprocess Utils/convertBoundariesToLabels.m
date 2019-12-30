function label_vol = convertBoundariesToLabels(boundary_points,vol_size,separate_bg)
% Create a label volume from boundary points - each layer is given a
% separate layer
%
% boundary points is an IxJxN array, where I and J are the x and y
% directions of the OCT data and N is the number of boundaries. Each layer
% is defined between consecutive boundaries so there are N-1 layers
% labeled.
%
% Background pixels (anything above or below the first and last boundary,
% respectively) are given a label of 0 unless separate_bg == true, in which
% case a separate label is used for background above and below the retina;
% above will have label 0 while below will have label N since N-1 is the
% number of labeled layers.

if nargin < 3
    separate_bg = false;
end

if ~isequal(vol_size(2:3),size(boundary_points(:,:,1)))
    error('input vol_size must match size of boundary points!')
end

% Construct label volume
% Boundary goes at bottom of layer
pts_r = round(boundary_points);
label_vol = zeros(vol_size,'uint8');
for ii = 1:size(label_vol,2)
    for jj = 1:size(label_vol,3)
        for kk = 1:(size(pts_r,3)-1)
            label_vol((pts_r(ii,jj,kk)+1):(pts_r(ii,jj,kk+1)),ii,jj) = kk;
        end
        if separate_bg
            label_vol((pts_r(ii,jj,end)+1):end,ii,jj) = size(pts_r,3);
        end
    end
end