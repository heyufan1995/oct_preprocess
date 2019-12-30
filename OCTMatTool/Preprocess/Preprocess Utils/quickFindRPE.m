function [rpe_bd, ilm_bd, shifts, retina_mask, bds, nbpt] = quickFindRPE(img_vol,voxel_size)
% Quick boundary detection by downsampling (doesn't have to be exact)

% Andrew Lang
% $Id: quickFindRPE.m,v 1.1 2015/01/12 18:52:29 andrew Exp $

if isstruct(voxel_size)
    % Header
    voxel_size = [voxel_size.ScaleZ voxel_size.ScaleX voxel_size.Distance];
end

vol_size = size(img_vol);
if numel(vol_size) == 2
    vol_size(3) = 1;
end

% Down sample to target resolution using nearest neighbor
voxel_size_target = [0.0039 0.0222 0.2443];
voxel_size_target = [0.0039 0.0222*3 voxel_size(3)];  
voxel_scale = round(voxel_size_target./voxel_size);

new_size = vol_size;
new_size = round(new_size./voxel_scale);
img_vol = imresize(img_vol,new_size(1:2));
voxel_size = voxel_scale.*voxel_size;

% Detect retina boundaries
header.ScaleX = voxel_size(2);
header.ScaleZ = voxel_size(1);
header.Distance = voxel_size(3);
[retina_mask,shifts,bds,nbpt] = retinaDetector(img_vol,header,'cirrus_sm',false);    

shifts = voxel_scale(1)*imresize(shifts,[vol_size(2) vol_size(3)]);
if nargout > 3
    bds = voxel_scale(1)*imresize(bds,[vol_size(2) vol_size(3)]);
    retina_mask = imresize(retina_mask,[vol_size(1) vol_size(2)],'nearest');
end

% RPE boundary resized to original data size
rpe_bd = voxel_scale(1)*bds(:,:,3);
rpe_bd = imresize(rpe_bd,[vol_size(2) vol_size(3)]);
rpe_bd(rpe_bd < 1) = 1;
rpe_bd(rpe_bd > vol_size(1)) = vol_size(1);

ilm_bd = voxel_scale(1)*bds(:,:,1);
ilm_bd = imresize(ilm_bd,[vol_size(2) vol_size(3)]);
ilm_bd(ilm_bd < 1) = 1;
ilm_bd(ilm_bd > vol_size(1)) = vol_size(1);