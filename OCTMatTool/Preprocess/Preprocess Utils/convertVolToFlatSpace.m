function [oct_vol_flat,flat_boundary_inds,reg_boundaries] = convertVolToFlatSpace(oct_vol,header,flat_params)
% method = 1 - flatten to average boundary positions
% method = 2 - flatten to regression learned boundary positions
%
% To add? Flatten to input boundary positions (flatten to each boundary
% given and linearly interpolate between boundaries)

% Flat space is simply a transformation of the depth position at each voxel
%   Instead of transforming the flat space data to native space directly,
%   we convert the native coordinates to flat space and interpolate using
%   the inverse mapping

interp_type = flat_params.interp_method;

%% Get size of flat space layers
% (which will give us the flat space boundary positions)

reg_boundaries = thicknessRegression(flat_params,header);

% Thickness of regression layers
N = getFlatLayerThicknesses(flat_params,reg_boundaries);

%%

% Position of each voxel (depth direction)
Nsum = sum(N) + size(reg_boundaries,3);
vol_size = double([Nsum size(oct_vol,2) size(oct_vol,3) size(oct_vol,4)]);
volPoints = repmat((1:vol_size(1))',[1 vol_size(2:3)]);

% Transform points from flat space
volPoints = permute(volPoints,[2 3 1]);

flat_params.pixels_per_layer = N;
volPointsTF = convertBoundariesFromFlatSpace(volPoints,header,flat_params,reg_boundaries);
volPointsTF = permute(volPointsTF,[3 1 2]);

% Do interpolation
cl = class(oct_vol);
oct_vol = double(oct_vol);

oct_vol_flat = zeros(vol_size);

% Grid points of input and output data
[Y,X] = ndgrid(1:size(oct_vol,1),1:size(oct_vol,2));
[~,Xi] = ndgrid(1:size(oct_vol_flat,1),1:size(oct_vol_flat,2));

if flat_params.verbose > 1
    fprintf('Interpolating data...')
end
for i = 1:size(oct_vol,4)  
    for j = 1:size(oct_vol,3)      
        % Interpolate volume at grid points    
        F = griddedInterpolant(Y,X,oct_vol(:,:,j,i),interp_type);
        F.ExtrapolationMethod = 'none';
        oct_vol_flat(:,:,j,i) = F(volPointsTF(:,:,j),Xi);
    end
    if flat_params.verbose > 1
        fprintf('%.0f%% ',i/size(oct_vol,4)*100);
    end
end
if flat_params.verbose > 1
    fprintf('\n')
end

oct_vol_flat = cast(oct_vol_flat,cl);

N_csum = cumsum(cat(1,0,N(1:end-1))) + (0:(length(N)-1))';
flat_boundary_inds = N_csum(2:end);