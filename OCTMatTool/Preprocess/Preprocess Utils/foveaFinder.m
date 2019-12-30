function [f_cen, th_val] = foveaFinder(bds,header,smoothfit,template_match,th_input,flip_template,legacy_sm)
% Compute the location of the center of the fovea computed generally as the
% point where the total retina thickness is smallest. Two methods can be
% used: the minimum thickness or template matching of a centered healthy
% control thickness template using cross correlation. If the thickness is
% used, the fovea point is found as that having the minimum thickness after
% averaging all thicknesses in a 300 um radius of each pixel .
%
% Inputs:
%   bds - an X x Z x S array of boundary points
%         an X x Z x (S-1) array of layer thicknesses
%         a Y x X x Z binary mask of the retina pixels
%
%         Y = # of pixels per A-scan
%         X = # of A-scans per B-scan, 
%         Z = # of B-scans
%         S = # of boundaries
%   header - OCT header with pixel sizes
%   smoothfit - flag to fit a cubic smoothing spline to the thickness map
%       used to estimate the fovea center (more robust?)
%   template_match - flag to use template matching to find the fovea center
%   th_input - flag stating whether or not layer thicknesses as opposed to
%       boundary points are input (ignored if input is a binary volume)
%   flip_template - flip the center point so that it matches the fundus
%       image view
%   legacy_sm - if smoothfit, estimate fovea using old method (for
%       consistency of thicknessStatistics function)
%
% Outputs:
%   f_cen - the center point of the fovea
%   th_val - retina thickness at this point
%
% Alternate methods to implement (?) include minimum thickness overall, and
% minimum thickness based on some sort of quadratic fitting

% Updates:
%   1/11/2016 - Use fminsearch to find minimum of spline surface fit when
%               using smoothfit = true (slightly faster, more accurate)

% Andrew Lang

if nargin < 3
    smoothfit = false;
end
if nargin < 4
    template_match = false;
end
if nargin < 5
    th_input = false;
end
if nargin < 6
    flip_template = false;
end
if nargin < 7
    % Need to test this more before moving forward with the new method
    legacy_sm = true;
end

if ~isfield(header,'angle')
    header.angle = 0;
end

scaleX = header.ScaleX*1000;
scaleY = header.Distance*1000;
scaleZ = header.ScaleZ*1000;

% Compute retina thickness
if islogical(bds)
    % Binary volume
    th = squeeze(sum(double(bds),1))*scaleZ;
else
    % Look at total distance between top and bottom boundary (probably ILM
    % and BM boundaries)
    if ~th_input
        th = sum(diff(bds,1,3),3)*scaleZ;
    else
        % Inputs are thickness maps
        th = sum(bds,3)*scaleZ;
    end
end

if isvector(th)
    if size(th,2) > size(th,1)
        th = th';
    end
    
    % Single b-scan, need to implement these methods in 2D if we want to
    % use them...
    smoothfit = false;
    template_match = false;
end

if flip_template
    th = fliplr(th');
end

th_orig = th;
if template_match
    % Template matching using retinal thickness template within 0.75 mm radius of fovea
    
    if any(isnan(th(:)));
        th = inpaint_nans(th);
    end
    
    % Get average retina thickness map
    S = load('mean_th');
    
    % Flip if right eye
    if ~strncmp(header.ScanPosition,'OD',2)
        S.mean_th = flipdim(S.mean_th,1);
    end
    mean_th = S.mean_th;
    % Extract center region around fovea
    i_cen = round(size(S.mean_th)/2);
    % Rotate if necessary
    if abs(header.angle) > 0.01
        x = ((1:size(mean_th,2))-i_cen(2))*S.scale_mean(3);
        y = ((1:size(mean_th,1))-i_cen(1))*S.scale_mean(2);
        [x,y] = meshgrid(x,y);
        [xr,yr] = rot2d(x,y,header.angle);
        mean_th = interp2(x,y,mean_th,xr,yr);
    end    
    r = 0.75;
    ybds = (i_cen(1)-round(r/S.scale_mean(2))):(i_cen(1)+round(r/S.scale_mean(2)));
    xbds = (i_cen(2)-round(r/S.scale_mean(3))):(i_cen(2)+round(r/S.scale_mean(3)));
    mean_th = mean_th(ybds,xbds);
    % Rescale to data resolution
    scale_ratio = S.scale_mean(2:3)./[header.ScaleX header.Distance];
    new_size = round(scale_ratio.*((size(mean_th)-1)/2))*2+1;
    mean_th = imresize(mean_th,new_size);
    % Normalize cross correlation matching
    xc = -circshift(normxcorr2(mean_th,th),(size(mean_th)-1)/2);
    xc(1:(size(mean_th,1)-1),:) = []; % crop to size
    xc(:,1:(size(mean_th,2)-1)) = [];

    th = xc;    
else
    % Get average thickness in 300 um radius circle
    r = 300;
    nx = ceil(r/scaleX); ny = ceil(r/scaleY);
    py = (-ny:ny)*scaleY; px = (-nx:nx)*scaleX;
    [Y,X] = meshgrid(py,px);
    R = sqrt(X.^2+Y.^2);
    m = R < r;
    th = imfilter(th,double(m)/sum(m(:)),'replicate');
end

% Only check the center part of the volume (assume square-ish volume)
cr1 = round([0.35 0.65]*size(th,2));
cr2 = round([0.35 0.65]*size(th,1));
cr1(cr1 < 1) = 1; % single b-scan
th_cen = th(cr2(1):cr2(2),cr1(1):cr1(2));

if smoothfit
    % Fit a smoothing spline and find the minimum that way
    x = {(1:size(th_cen,1))*scaleX (1:size(th_cen,2))*scaleY};
    p = csaps(x,th_cen);
    
    if legacy_sm
        % Minimum by dense sampling of spline (1000x1000 fundus map)
        x = {linspace(scaleX,size(th_cen,1)*scaleX,1000),linspace(scaleY,size(th_cen,2)*scaleY,1000)};
        sp_fit = fnval(p,x);    
        [~, ptm] = min(sp_fit(:));
        [ptm1, ptm2] = ind2sub(size(sp_fit),ptm);

        % Scale back to pixel space
        ptm1 = x{1}(ptm1)/scaleX;
        ptm2 = x{2}(ptm2)/scaleY;
    else
        % Find minimum of fitted spline
        x2 = fminsearch(@(x) fnval(p,{x(1),x(2)}),...
                        [size(th_cen,1)/2*scaleX size(th_cen,2)/2*scaleY]);

        % Scale back to pixel space
        ptm1 = x2(1)/scaleX;
        ptm2 = x2(2)/scaleY;    
    end
else
    % Fovea center found from the minumum
    [~, ptm] = min(th_cen(:));
    [ptm1, ptm2] = ind2sub(size(th_cen),ptm);
end

% Fovea center adding back the crop
f_cen = [ptm1 + cr2(1) - 1, ptm2 + cr1(1) - 1];
th_val = th_orig(round(f_cen(1)),round(f_cen(2)));

if flip_template
    f_cen(2) = size(th_orig,2)-f_cen(2)+1;
    f_cen = fliplr(f_cen);
end
