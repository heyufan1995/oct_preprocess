function [int_vals_sm, sp_fit, bscan_fit] = intensityCorrectionFactor(img_vol,bdry,voxel_size,minv,maxv,above,maxint)
% Compute the intensity correction factor used by the normalizeOCTVolume
% function

% Andrew Lang
% $Id: intensityCorrectionFactor.m,v 1.1 2015/01/12 18:52:25 andrew Exp $

if above
    bdry(bdry-maxv < 1) = maxv + 1;
else
    bdry(bdry+maxv > size(img_vol,1)) = size(img_vol,1) - maxv;
end

ef_mask = false(size(img_vol));
for i = 1:size(img_vol,2)
    for j = 1:size(img_vol,3)
        if above
            ef_mask((bdry(i,j,1)-maxv):(bdry(i,j,1)-minv),i,j) = true;
        else
            ef_mask((bdry(i,j,1)+minv):(bdry(i,j,1)+maxv),i,j) = true;
        end
    end
end

if maxint
    img_vol_masked = img_vol;
    img_vol_masked(~ef_mask) = 0;
    int_vals = squeeze(max(img_vol_masked,[],1));
else
    img_vol_masked = 1-img_vol;
    img_vol_masked(img_vol_masked == 1) = nan; % for zero background values (should be at least some noise)
    img_vol_masked(~ef_mask) = nan;
    int_vals = squeeze(nanmedian(img_vol_masked,1));
end

% Fill in nans if there are any
int_vals_nv = int_vals;
if any(isnan(int_vals(:)))
    xpts = 1:size(int_vals,1);
    for i = 1:size(int_vals,2)
        nv = isnan(int_vals(:,i));
        if any(nv)
            if sum(nv) < 0.5*size(int_vals,1)
                nvpts = median(int_vals(~nv,i));
%                     nvpts = interp1(xpts(~nv),ilm_int(~nv,i),xpts(nv),'nearest','extrap');
            else
                nvpts = nanmedian(int_vals(:));
            end
            int_vals_nv(nv,i) = nvpts;
        end
    end
end

% Start with a smooth surface fit
x = {(1:size(int_vals,1))*voxel_size(2) (1:size(int_vals,2))*voxel_size(3)};
[sp_fit,p] = csaps(x,double(int_vals_nv),0.9,x);

int_vals2 = int_vals./sp_fit;
int_vals2_nv = int_vals_nv./sp_fit;

% Now fit a line to the intensity profile of each A-scan
s = warning;
warning('off');
bscan_fit = zeros(size(int_vals));
for i = 1:size(img_vol,3)
    x = (1:size(int_vals,1))'*voxel_size(2);
%             X = [x x.^2 x.^3 x.^4];
    X = x;

    if sum(isnan(int_vals2(:,i))) > 0.5*size(int_vals,1)
        st = robustfit(X,int_vals2_nv(:,i),'bisquare',2);
    else
        st = robustfit(X,int_vals2(:,i),'bisquare',2);
    end
    bscan_fit(:,i) = polyval(flipud(st),(1:size(img_vol,2))*voxel_size(2));

%         figure,plot((1:size(img_vol,2))*voxel_size(2),int_vals2(:,i),'r')
%         hold on,plot((1:size(img_vol,2))*voxel_size(2),int_vals_sm(:,i),'b')
end
warning(s);

if maxint
    int_vals_sm = bscan_fit.*sp_fit;
else
    int_vals_sm = 1-bscan_fit.*sp_fit;
end
