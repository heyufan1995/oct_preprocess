function slo = octFundus(img_vol,rpe_bds,noflip,header)
% Generate an OCT fundus projection image
%   If no RPE boundary is input then the output is just the sum of
%   intensities along each A-scan, otherwise, it is the sum of hte
%   intensities from the RPE to 20 pixels above it.
%
% Added support for no flip - set 'noflip' to true if you don't want the
% output to be in the fundus orientation

% Andrew Lang
% $Id: octFundus.m,v 1.1 2015/01/12 18:54:06 andrew Exp $

if nargin < 3
    noflip = false;
end
if nargin < 4
    header.ScaleZ = 3.8/1000;
end

if ~isempty(rpe_bds)
    % Average 80 microns above the RPE boundary
    ef_th = round(80/(header.ScaleZ*1000));
    rpe_bds = round(rpe_bds);
    rpe_bds_top = rpe_bds - ef_th;
    rpe_bds_top(rpe_bds_top < 1) = 1;
    ef_mask = false(size(img_vol));
    for i = 1:size(img_vol,2)
        for j = 1:size(img_vol,3)
            ef_mask(rpe_bds_top(i,j,1):rpe_bds(i,j,1),i,j) = true;
        end
    end
    img_vol(~ef_mask) = 0;
    if noflip
        slo = squeeze(sum(img_vol,1));
        
        % Normalize each B-scan then across B-scans 
        slo = bsxfun(@rdivide,slo,mean(slo,1));
        slo = bsxfun(@rdivide,slo,mean(slo,2));
    else
        slo = flipud(squeeze(sum(img_vol,1))');
        
        % Normalize each B-scan then across B-scans 
        slo = bsxfun(@rdivide,slo,mean(slo,2));
        slo = bsxfun(@rdivide,slo,mean(slo,1));
    end    
    
    % Adjust contrast
    slo = (slo-0.5)/(1.2-0.5);
    slo(slo < 0) = 0;
    slo(slo > 1) = 1;
else
    % Don't have RPE boundary so just average along A-scans
    if noflip
        slo = squeeze(sum(img_vol,1));
        
        % Normalize each B-scan then across B-scans 
        slo = bsxfun(@rdivide,slo,mean(slo,1));
        slo = bsxfun(@rdivide,slo,mean(slo,2));
    else
        slo = flipud(squeeze(sum(img_vol,1))');
        
        % Normalize each B-scan then across B-scans 
        slo = bsxfun(@rdivide,slo,mean(slo,2));
        slo = bsxfun(@rdivide,slo,mean(slo,1));
    end
end