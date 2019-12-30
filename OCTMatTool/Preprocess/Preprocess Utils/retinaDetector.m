function [retinaMask, shifts, boundaries, nbpt] = retinaDetector(img_vol,header,paramSet,doplots)
% Detect the retina boundaries. The retina mask will have a value of 0
% outside of the retina, a value of  1 between the ilm and inner segment,
% and a value of 2 between the inner segment and Bruch's membrane (bottom
% of the RPE)
%
% Note, code was not optimized to accurately find the ISOS boundary, only
% the ILM and BrM
%
% Shifts is an YxZ vector containing the number of pixels to shift each
% a-scan to flatten the image. Use retinaFlatten.m to flatten the data
% using these shifts.
%
% Boundaries is an YxZx3 array containing the value of the ILM, ISOS and
% BM boundaries

% Fixes: %   - zero values on edges caused a problem
%   - enforce layer ordering
%   - boundaries out of the volume
%   - cirrus b-scan independent processing
%   - improved BM detection by using depth weighting
%   - changed filter boundary conditions from 'symmetric' to 'replicate'
%   - BM now fit with quadratic for bsc_indep setting
%   - More conservative outlier detection makes for larger nbpt counts

% File formerly named retinaDetector2_scale.m

% Andrew Lang
% $Id: retinaDetector.m,v 1.2 2015/01/12 18:52:25 andrew Exp $

if nargin < 2
    error('number of arguments must be at least 2!')
elseif nargin < 3
    paramSet = 'default';
    doplots = false;
elseif nargin < 4
    doplots = false;
end

params = getParams(paramSet);

sl = round(size(img_vol,3)/2);

% Maximum distance from ILM to ISOS
maxdist = params.maxdist;
% Maximum distance from ISOS to BM
maxdist_bm = 116.02; % ~30 pixels in spectralis
% Minimum distance from ISOS to BM
isosThresh = 20; 
% Median filter outlier threshold distance and kernel
dc_thresh = 10;
mf_k = 140;

% Process B-scans independently
bsc_indep = params.bsc_indep;
%% vertical scans! added Bhavna J Antony, 08/26/16
%% since vertical scans were showing more motion artifacts
%% than horizontal scans. Can be removed if/when Spectralis
%% fixes the issue 
if isfield(header, 'angle')
    if abs(abs(header.angle)-90) < 25
        bsc_indep = 1;
    end
end

% Sigma values for smoothing final surfaces
% Through plane direction
sigma_tp_ilm = 91.62; % ~0.75 pixels in spectralis
sigma_tp_isos = 91.62; % ~0.75 pixels in spectralis
sigma_tp_bm = 244.32; % ~2 pixels in spectralis
% Lateral direction
sigma_lat_ilm = params.sigma_lat_ilm; % ~10 pixels in spectralis
sigma_lat_isos = params.sigma_lat_isos; % ~10 pixels in spectralis
sigma_lat_bm = params.sigma_lat_bm; % ~20 pixels in spectralis

%% Convert all values from micron to pixel
sigma_lat = params.sigma_lat/(header.ScaleX*1000);
sigma_ax = params.sigma_ax/(header.ScaleZ*1000);
distConst = round(params.distconst/(header.ScaleZ*1000));
maxdist = round(maxdist/(header.ScaleZ*1000));
maxdist_bm = round(maxdist_bm/(header.ScaleZ*1000));
isosThresh = round(isosThresh/(header.ScaleZ*1000));
dc_thresh = round(dc_thresh/(header.ScaleZ*1000)*(128/6)*header.Distance); % scale by b-scans distance
% mf_k = round(mf_k/(header.Distance*1000))*2+1;
mf_k = [round(mf_k/(header.ScaleX*1000)) round(mf_k/(header.Distance*1000))]*2+1;
sigma_tp_ilm = sigma_tp_ilm/(header.Distance*1000);
sigma_tp_isos = sigma_tp_isos/(header.Distance*1000);
sigma_tp_bm = sigma_tp_bm/(header.Distance*1000);
sigma_lat_ilm = sigma_lat_ilm/(header.ScaleX*1000);
sigma_lat_isos = sigma_lat_isos/(header.ScaleX*1000);
sigma_lat_bm = sigma_lat_bm/(header.ScaleX*1000);

%% Handle zero or nan nan values on the borders

% The spectralis often has zero or nan values around the edges of the image
% due the eye tracking. Fill in these values with the nearest pixel

img_vol(isnan(img_vol)) = 0;

% Fill in from the left side
[~,inds] = max(img_vol>0,[],2); % Index of first nonzero value
for i = 1:size(img_vol,1)
    for j = 1:size(img_vol,3)
        p = inds(i,1,j);
        if p > 1 && p < i
            if p < (size(img_vol,2)-2)
                % Avoid using low intensity edge pixels
                img_vol(i,1:(p+1),j) = img_vol(i,p+2,j); 
            else
                img_vol(i,1:(p-1),j) = img_vol(i,p,j);
            end
        end
    end
end
% Fill in from the right side
[~,inds] = max(flipdim(img_vol>0,2),[],2); % Index of last nonzero value
inds = size(img_vol,2) - inds + 1;
for i = 1:size(img_vol,1)
    for j = 1:size(img_vol,3)
        p = inds(i,1,j);
        if p < size(img_vol,2) && size(img_vol,2)-p < i
            if p > 2
                % Avoid using low intensity edge pixels
                img_vol(i,(p-1):end,j) = img_vol(i,p-2,j); 
            else
                img_vol(i,(p+1):end,j) = img_vol(i,p,j);
            end
        end
    end
end
% Fill in from the top
mv = mean(img_vol(:));
[~,inds] = max(img_vol>0,[],1); % Index of first nonzero value
for i = 1:size(img_vol,2)
    for j = 1:size(img_vol,3)
        p = inds(1,i,j);
        if p > 1
            if p < (size(img_vol,1)-2)
                % Avoid using low intensity edge pixels
                if img_vol(p+2,i,j) < mv
                    img_vol(1:(p+1),i,j) = img_vol(p+2,i,j);
                else
                    % Probably cut through the retina so keep a gradient
                    % here
                    img_vol(1:(p+1),i,j) = mv;
                end
            else
                img_vol(1:(p-1),i,j) = img_vol(p,i,j);
            end
        end
    end
end
% Fill in from the bottom
[~,inds] = max(flipdim(img_vol>0,1),[],1); % Index of first nonzero value
inds = size(img_vol,1) - inds + 1;
for i = 1:size(img_vol,2)
    for j = 1:size(img_vol,3)
        p = inds(1,i,j);
        if p < size(img_vol,1)
            if p > 2
                % Avoid using low intensity edge pixels
                img_vol((p-1):end,i,j) = img_vol(p-2,i,j); 
            else
                img_vol((p+1):end,i,j) = img_vol(p,i,j);
            end
        end
    end
end

if doplots
    h1 = figure;
    imagesc(img_vol(:,:,sl)), colormap gray
end

%% Pre-processing

% Gaussian filter data
% grad = calculateFeatures2D({'anigauss',{sigma_ax,sigma_lat,0,0,[]}},img_vol);
grad = imfilter(img_vol,fspecial('gaussian', [2*round(2*sigma_ax)+1 1], sigma_ax),'replicate');
grad = imfilter(grad,fspecial('gaussian', [1 2*round(2*sigma_lat)+1], sigma_lat),'replicate');

if doplots
    figure,imagesc(grad(:,:,sl)), colormap gray
end

% Gradient of the image in the y-direction
grad = -imfilter(grad,fspecial('sobel'),'replicate');

if doplots
    figure,imagesc(grad(:,:,sl)), colormap gray
end

%% Find ILM and ISOS boundaries

grad_o = grad;

% Largest gradient, this is either the top layer or the center layer
[~, max1pos] = max(grad,[],1);
max1pos = squeeze(max1pos);
if isvector(max1pos)
    max1pos = max1pos';
end

% Find the second largest gradient to the max gradient at a distance of at
% least distConst away but not more than maxdist away
for i = 1:size(grad,2)
    for j = 1:size(grad,3)
        % min distance
        dc = distConst;
        if (max1pos(i,j)-distConst) < 1
            dc = max1pos(i,j) - 1;
        elseif (max1pos(i,j)+distConst) > size(grad,1)
            dc = size(grad,1) - max1pos(i,j);
        end        
        grad((max1pos(i,j)-dc):(max1pos(i,j)+dc),i,j) = 0;
        
        % max distance
        if (max1pos(i,j)-maxdist) > 0
            grad(1:(max1pos(i,j)-maxdist),i,j) = 0;
        end
        if (max1pos(i,j)+maxdist) <= size(grad,1)
            grad((max1pos(i,j)+maxdist):end,i,j) = 0;
        end
    end
end

[~,max2pos] = max(grad,[],1);
max2pos = squeeze(max2pos);
if isvector(max2pos)
    max2pos = max2pos';
end

ilm = min(max1pos,max2pos);
isos = max(max1pos,max2pos);

if doplots
    figure(h1),hold on,plot(1:size(img_vol,2),ilm(:,sl),'r.','markersize',5)
    figure(h1),hold on,plot(1:size(img_vol,2),isos(:,sl),'b.','markersize',5)
end

%% Find BM boundary

grad = grad_o;

% BM is largest negative gradient below the ISOS
for i = 1:size(grad,2)
    for j = 1:size(grad,3)  
        grad(1:isos(i,j)+isosThresh,i,j) = 0;
        
        if (isos(i,j)+maxdist_bm) <= size(grad,1)
            grad((isos(i,j)+maxdist_bm):end,i,j) = 0;
        end
    end
end

% To encourage boundary points closer to the top of the image, weight
% linearly by depth (sometimes the BM boundary gets confused with random
% gradients in the choroid, so this will alleviate that problem)
lin = bsxfun(@plus,(1:size(grad,1))',size(grad,1)-(shiftdim(isos,-1)+maxdist_bm));
lin = -0.5/size(grad,1)*lin+1;
grad = grad.*lin;

[~,bot] = min(grad,[],1); 
bm = squeeze(bot);

if isvector(bm)
    bm = bm';
end

if doplots
    figure(h1),hold on,plot(1:size(img_vol,2),bm(:,sl),'g.','markersize',5)
end

%% Detect outliers

if bsc_indep
    % Detect outliers from total retina thickness
    th = bm-ilm;
%     th_med = medfilt2(th,[1 mf_k],'symmetric');
    th_med = medfilt2(th,mf_k,'symmetric');
    
    bpt = (abs(th-th_med) > dc_thresh);
    
%     % Also remove points with large gradient in the B-scan direction as a
%     % little added protection against random spikes (a random spike could
%     % have a thickness that looks reasonable)
%     ilm_grad = abs(gradient(ilm')') > 3*dc_thresh;
%     bm_grad = abs(gradient(bm')') > 3*dc_thresh;
%     
%     % Dilate by 1 pixel in B-scan direction since the gradient of the spike
%     % is 0 at the center (gradient([0 1 0] = [1 0 -1])
%     ilm_grad = imdilate(ilm_grad,strel('line',3,90));
%     bm_grad = imdilate(bm_grad,strel('line',3,90));
%     
%     bpt = (abs(th-th_med) > dc_thresh) | ilm_grad | bm_grad;
else
    % Median filter surfaces to detect outliers
%     ilm_med = medfilt2(ilm,[1 mf_k],'symmetric');
%     isos_med = medfilt2(isos,[1 mf_k],'symmetric');
%     bm_med = medfilt2(bm,[1 mf_k],'symmetric');
    
    ilm_med = medfilt2(ilm,mf_k,'symmetric');
    isos_med = medfilt2(isos,mf_k,'symmetric');
    bm_med = medfilt2(bm,mf_k,'symmetric');
    
    bpt = (abs(ilm-ilm_med) > dc_thresh) | (abs(isos-isos_med) > dc_thresh) ...
          | (abs(bm-bm_med) > dc_thresh);
end


% Fill in outlier points
ilm(bpt) = nan;
isos(bpt) = nan;
bm(bpt) = nan;
nbpt = 0;
if any(bpt(:))
    nbpt = sum(bpt(:));
    if bsc_indep
        % Fit a quadratic to BM
        x = (1:size(ilm,1))';
        for j = 1:size(ilm,2)            
            p = polyfit(x(~bpt(:,j)),bm(~bpt(:,j),j),3);
            bm(:,j) = polyval(p,x);
        end
        
        % Linearly interpolate ILM and ISOS
        nv = any(isnan(ilm));
        xpts = 1:size(ilm,1);
        for j = 1:size(ilm,2)
            if nv(j)
                nv2 = ~isnan(ilm(:,j));
                ilm(:,j) = interp1(xpts(nv2),ilm(nv2,j),xpts,'linear','extrap');
                isos(:,j) = interp1(xpts(nv2),isos(nv2,j),xpts,'linear','extrap');
%                 bm(:,j) = interp1(xpts(nv2),bm(nv2,j),xpts,'linear','extrap');
            end
        end
    else
        ilm = inpaint_nans(ilm);
        isos = inpaint_nans(isos);
        bm = inpaint_nans(bm);
    end
end

%% Get final boundaries by smoothing

% Finally smooth surfaces
if ~bsc_indep
    ilm = imfilter(ilm,fspecial('gaussian',[1 2*round(3*sigma_tp_ilm)+1],...
                                sigma_tp_ilm),'replicate');
    isos = imfilter(isos,fspecial('gaussian',[1 2*round(3*sigma_tp_isos)+1],...
                                  sigma_tp_isos),'replicate');
    bm = imfilter(bm,fspecial('gaussian',[1 2*round(3*sigma_tp_bm)+1],...
                              sigma_tp_bm),'replicate');
    bm = imfilter(bm,fspecial('gaussian',[2*round(3*sigma_lat_bm)+1 1],...
                          sigma_lat_bm),'replicate');     
end

ilm = imfilter(ilm,fspecial('gaussian',[2*round(3*sigma_lat_ilm)+1 1],...
                            sigma_lat_ilm),'replicate');
isos = imfilter(isos,fspecial('gaussian',[2*round(3*sigma_lat_isos)+1 1],...
                              sigma_lat_isos),'replicate');
% bm = imfilter(bm,fspecial('gaussian',[2*round(3*sigma_lat_bm)+1 1],...
%                           sigma_lat_bm),'replicate');
 
% Enforce ordering and a very small minimum thickness
bmilm = (bm-ilm)*header.ScaleZ*1000 < 100;
ilm(bmilm) = bm(bmilm) - 100/header.ScaleZ/1000;
bmisos = (bm-isos)*header.ScaleZ*1000 < 10;
isos(bmisos) = bm(bmisos) - 10/header.ScaleZ/1000;
isosilm = (isos-ilm)*header.ScaleZ*1000 < 90;
isos(isosilm) = ilm(isosilm) + 90/header.ScaleZ/1000;

if doplots
    figure(h1),hold on,plot(1:size(img_vol,2),ilm(:,sl),'m.','markersize',5)
    figure(h1),hold on,plot(1:size(img_vol,2),isos(:,sl),'c.','markersize',5)
    figure(h1),hold on,plot(1:size(img_vol,2),bm(:,sl),'y.','markersize',5)
end

% Make sure we are not out of the volume
ilm(ilm<1) = 1;
ilm(ilm>size(img_vol,1)) = size(img_vol,1);
isos(isos<1) = 1;
isos(isos>size(img_vol,1)) = size(img_vol,1);
bm(bm<1) = 1;
bm(bm>size(img_vol,1)) = size(img_vol,1);

% Create mask volume
retinaMask = zeros(size(img_vol),'uint8');
for i = 1:size(img_vol,2)
    for j = 1:size(grad,3)
        retinaMask(round(ilm(i,j)):round(isos(i,j)-1),i,j) = 1;
        retinaMask(round(isos(i,j)):round(bm(i,j)),i,j) = 2;
    end
end

boundaries = cat(3,ilm,isos,bm);

% Shifts (we shift so that the bottom layer is in the center of the image)
% shifts = bm' - mean(bm) - (round(size(image,1)/2)-mean(bm));
shifts = bsxfun(@minus,bm,mean(bm,1)+(round(size(img_vol,1)/2)-mean(bm,1)));


function params = getParams(paramSet)

switch lower(paramSet)
    case {'default','spectralis','hc','mme'}
        % Originally for the spectralis
        params.sigma_lat = 16.67;
        params.sigma_ax = 11.6;
        params.distconst = 96.68;
        params.sigma_lat_ilm = 55.56;
        params.sigma_lat_isos = 55.56;
        params.sigma_lat_bm = 111.13;
        params.maxdist = 386.73; % ~100 pixels in spectralis
        params.bsc_indep = false;
     case {'dme'}
        % Originally for the spectralis
        params.sigma_lat = 16.67;
        params.sigma_ax = 11.6;
        params.distconst = 96.68;
        params.sigma_lat_ilm = 55.56;
        params.sigma_lat_isos = 55.56;
        params.sigma_lat_bm = 111.13;
        params.maxdist = 386.73*2; % ~100 pixels in spectralis
        params.bsc_indep = false;
    case 'cirrus'
        % Slightly different parameters for cirrus
        params.sigma_lat = 2*16.67;
        params.sigma_ax = 0.5*11.6;
        params.distconst = 96.68;
        params.sigma_lat_ilm = 55.56;
        params.sigma_lat_isos = 55.56;
        params.sigma_lat_bm = 111.13;
        params.maxdist = 386.73;
        params.bsc_indep = true;
    case 'cirrus_sm'
        % Additional smoothing for RPE
        params.sigma_lat = 2*16.67;
        params.sigma_ax = 0.5*11.6;
        params.distconst = 96.68;
        params.sigma_lat_ilm = 55.56;
        params.sigma_lat_isos = 55.56;
        params.sigma_lat_bm = 200;
        params.maxdist = 386.73;
        params.bsc_indep = true;
    case 'rp'
        % Retinitis pigmentosa data requires extra smoothing to remove
        % outlier bumps
        params.sigma_lat = 2*16.67;
        params.sigma_ax = 0.5*11.6;
        params.distconst = 50;
        params.sigma_lat_ilm = 200;
        params.sigma_lat_isos = 300;
        params.sigma_lat_bm = 200;
        params.maxdist = 386.73;
        params.bsc_indep = true;
    case 'phantom'
        params.sigma_lat = 5;
        params.sigma_ax = 5;
        params.distconst = 150;
        params.sigma_lat_ilm = 55.56;
        params.sigma_lat_isos = 55.56;
        params.sigma_lat_bm = 111.13;
        params.maxdist = 550;
        params.bsc_indep = false;
end