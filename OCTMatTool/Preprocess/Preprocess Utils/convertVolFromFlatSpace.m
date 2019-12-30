function oct_vol = convertVolFromFlatSpace(oct_vol_flat,header,flat_params)
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

retina_boundaries = flat_params.retina_boundaries;

% Position of each voxel (depth direction)
vol_size = double([header.SizeZ size(oct_vol_flat,2) size(oct_vol_flat,3) size(oct_vol_flat,4)]);
volPoints = repmat((1:vol_size(1))',[1 vol_size(2:3)]);

% Crop for efficiency
maxPt = round(max(retina_boundaries(:)))+ceil(flat_params.pad_below/(header.ScaleZ*1000));
minPt = round(min(retina_boundaries(:)))-ceil(flat_params.pad_above/(header.ScaleZ*1000));
if maxPt > vol_size(1)
    maxPt = vol_size(1);
end
if minPt < 1
    minPt = 1;
end
volPointsCrop = volPoints(minPt:maxPt,:,:);
crop_size = [size(volPointsCrop) size(oct_vol_flat,4)];
clear volPoints

% Transform points to flat space
volPointsCrop = permute(volPointsCrop,[2 3 1]);
volPointsTF = convertBoundariesToFlatSpace(volPointsCrop,header,flat_params);
volPointsTF = permute(volPointsTF,[3 1 2]);

% Do interpolation
cl = class(oct_vol_flat);
oct_vol_flat = double(oct_vol_flat);

oct_vol_crop = zeros(crop_size);

% Grid points of input and output data
[Y,X] = ndgrid(1:size(oct_vol_flat,1),1:size(oct_vol_flat,2));
[~,Xi] = ndgrid(1:crop_size(1),1:crop_size(2));

if flat_params.verbose > 1
    fprintf('Interpolating data...')
end
for i = 1:size(oct_vol_flat,4)
    for j = 1:size(oct_vol_flat,3)
        % Interpolate volume at grid points    
        F = griddedInterpolant(Y,X,oct_vol_flat(:,:,j,i),interp_type);
        F.ExtrapolationMethod = 'none';
        oct_vol_crop(:,:,j,i) = F(volPointsTF(:,:,j),Xi);
    end
    if flat_params.verbose > 1
        fprintf('%.0f%% ',i/size(oct_vol_crop,4)*100);
    end
end
if flat_params.verbose > 1
    fprintf('\n')
end

% Uncrop
oct_vol = nan(vol_size);
oct_vol(minPt:maxPt,:,:,:) = oct_vol_crop;

oct_vol = cast(oct_vol,cl);