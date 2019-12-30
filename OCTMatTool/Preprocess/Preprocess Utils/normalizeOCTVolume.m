function out_vol = normalizeOCTVolume(in_vol,method,header,retina_bds,downsamp)
% Different intensity normalization methods (method input)
%   1 - median filter standard deviation normalization
%   2 - median filter contrast stretching normalization (used in BOE paper)
%   3 - quantile contrast stretching
%   4 - Attenuation correction (Method of Girard et al. (2011))
%   5 - median filter contrast stretching with lower threshold
%   6 - normalize to RPE intensity
%   7 - normalize to RPE and vitreous intensities
%   8 - normalize to RPE and RNFL intensities
%   9 - normalize to smoothly varying bias field based on RPE intensity
%   10 - normalize to RPE intensity only on B-scans
%   11 - normalize using N3O
% Updates:
% Added option to input retina boundaries instead of computing them by
%   downsampling
% Reordered methods so they make more sense in order

% Andrew Lang
% $Id: normalizeOCTVolume.m,v 1.3 2015/01/12 18:52:25 andrew Exp $

if nargin < 5
    downsamp = false;
end
if nargin < 4
    retina_bds = [];
end
if nargin < 3
    header = [];
end
if nargin < 2
    method = 2;
end

if downsamp
    % Downsample data to fixed resolution
    % Note - does not affect all methods
    rs_res = [4 10];
    old_res = 1000*[header.ScaleZ header.ScaleX];
    rs_res = max(rs_res,old_res); % make sure we don't upsample
    
    sc = old_res./rs_res;
    new_size = round(size(in_vol(:,:,1)).*sc);
    in_vol_rs = imresize(in_vol,new_size);
else
    in_vol_rs = in_vol;
    rs_res = 1000*[header.ScaleZ header.ScaleX];
end

% Median filter kernel
hsz = 27;
h = [2*round(hsz/rs_res(1))+1 1];

out_vol = zeros(size(in_vol),class(in_vol));
if method == 1
    % Median filter standard deviation normalization
    
    % Get standard deviation of whole volume, we want each slice to have
    % approximately this standard deviation (which also keeps the
    % intensity range of 0-1 meaningful)
    stdv = nanstd(in_vol_rs(:));
    
    % Loop through each slice, median filter the image and normalize the
    % slice by the standard deviation of the median filtered image
    for j = 1:size(in_vol,3)
        med = medfilt2(in_vol_rs(:,:,j),h);
        out_vol(:,:,j) = in_vol(:,:,j)./nanstd(med(:))*stdv;
    end
elseif method == 2
    % Median filter contrast stretching normalization
    if isinteger(in_vol)
        mv = double(intmax(class(in_vol))); % imadjust requires double
    else
        mv = 1;
    end
    
    maxOffset = 0.05*mv;
    
    % Loop through each slice, median filter the image and normalize by
    % contrast stretching to the max of the median filtered image
    for j = 1:size(in_vol,3)
        med = medfilt2(in_vol_rs(:,:,j),h);
        ms = double(max(med(:)))+maxOffset;
        if ms > mv
            ms = mv;
        end
        out_vol(:,:,j) = imadjust(in_vol(:,:,j),[0 ms]/mv,[0 1]);
    end    
elseif method == 3
    % Quantile contrast stretching
    qle1 = 0.3;
    qle2 = 0.999;
    
    q1 = quantile(in_vol_rs(:),qle1);
    q2 = quantile(in_vol_rs(:),qle2);
    
    for j = 1:size(in_vol,3)
%         sl = in_vol(:,:,j);
%         q1 = quantile(sl(:),qle1);
%         q2 = quantile(sl(:),qle2);
        out_vol(:,:,j) = imadjust(in_vol(:,:,j),[q1 q2],[0 1]);
    end    
    
elseif method == 4
    % Method of Girard et al. (2011), Shadow removal and contrast
    % enhancement in OCT images of the human optic nerve head
    
    % Works on native data, not 4th root
    in_vol = in_vol.^2.^2; % faster than .^4
    
    cs = flipdim(cumtrapz(flipdim(in_vol,1),1),1);
    sigma = 15/(header.ScaleX*1000);
    cs = imfilter(cs,fspecial('gaussian',[1 2*round(3*sigma)+1],sigma),'symmetric');
    
    % Adaptive energy threshold - Mari et al (2014)
    E = flipdim(cumtrapz(flipdim(in_vol.^2,1),1),1);
    E = E < 0.0001;
    for i = 1:size(in_vol,3)
        ind = find(all(E(:,:,i),2),1,'first');
        cs(ind:end,:,i) = repmat(cs(ind,:,i),[size(cs,1)-ind+1 1]);
    end
    
    out_vol = in_vol./cs;
    out_vol(isnan(out_vol)) = min(out_vol(:)); % 0/0
    out_vol(isinf(out_vol)) = min(out_vol(:)); % 0/0
    
    out_vol = sqrt(sqrt(out_vol)); % faster than ^0.25;
    
elseif method == 5
    % Median filter contrast stretching normalization with lower threshold
    if isinteger(in_vol)
        mv = double(intmax(class(in_vol))); % imadjust requires double
    else
        mv = 1;
    end
    
    maxOffset = 0.05*mv;
    minOffset = 0.05;
    
    % Loop through each slice, median filter the image and normalize by
    % contrast stretching to the max of the median filtered image
    for j = 1:size(in_vol,3)
        med = medfilt2(in_vol_rs(:,:,j),h);
        ms = double(max(med(:)))+maxOffset;
        if ms > mv
            ms = mv;
        end
        mns = double(min(med(:)))+minOffset;
        out_vol(:,:,j) = imadjust(in_vol(:,:,j),[mns ms]/mv,[0 1]);
    end 

elseif method == 6
    % Normalize to RPE intensity

    voxel_size = [header.ScaleZ header.ScaleX header.Distance];
    
    if isempty(retina_bds)
        % Estimate RPE and ILM
        [rpe_bdry,~] = quickFindRPE(in_vol,voxel_size);
        rpe_bdry = round(rpe_bdry); 
    else
        rpe_bdry = round(retina_bds(:,:,2)); 
    end

    % Find the intensity correction factor for the RPE
    minv = 0;
    maxv = round(80/(voxel_size(1)*1000));    
    above = true;
    maxint = true;
    rpe_int_sm = intensityCorrectionFactor(in_vol,rpe_bdry,voxel_size,minv,maxv,above,maxint);
    
    % Protect against very small values
    rpe_int_sm(rpe_int_sm<0.05) = 0.05;    
    
    % Divide by correction factor to make RPE intensities approx. 1
    out_vol = bsxfun(@rdivide,in_vol,shiftdim(rpe_int_sm,-1))*.90;
    out_vol(out_vol>1) = 1;    
    
elseif method == 7
    % Normalize to RPE and vitreous intensities
    % Want to make RPE have a value of 1 and background a value of 0
                
    voxel_size = [header.ScaleZ header.ScaleX header.Distance];
    
    if isempty(retina_bds)
        % Estimate RPE and ILM
        [rpe_bdry,ilm_bdry] = quickFindRPE(in_vol,voxel_size);
        rpe_bdry = round(rpe_bdry); 
        ilm_bdry = round(ilm_bdry);
    else
        rpe_bdry = round(retina_bds(:,:,2)); 
        ilm_bdry = round(retina_bds(:,:,1));
    end
                
    % Find the intensity correction factor for the RPE
    minv = 0;
    maxv = round(80/(voxel_size(1)*1000));    
    above = true;
    maxint = true;
    rpe_int_sm = intensityCorrectionFactor(in_vol,rpe_bdry,voxel_size,minv,maxv,above,maxint);
    
    % Protect against very small values
    rpe_int_sm(rpe_int_sm<0.05) = 0.05; 
    
    % Divide by correction factor to make RPE intensities approx. 1
    out_vol = bsxfun(@rdivide,in_vol,shiftdim(rpe_int_sm,-1))*.90;
    out_vol(out_vol>1) = 1;

    % Now do a similar thing for the background intensities above the 
    %   ILM (divide 1-octvol by the fit of the background)
    minv = round(20/(voxel_size(1)*1000));
    maxv = round(80/(voxel_size(1)*1000));    
    above = true;
    maxint = false;
    ilm_int_sm = intensityCorrectionFactor(double(out_vol),ilm_bdry,voxel_size,minv,maxv,above,maxint);

    % Divide by correction factor to make background intensities approx. 1
    out_vol = 1-bsxfun(@rdivide,1-out_vol,shiftdim(1-ilm_int_sm,-1))*0.9;
%     out_vol = 1-bsxfun(@rdivide,1-out_vol,shiftdim(1-ilm_int_sm,-1));
    out_vol(out_vol<0) = 0;
    
elseif method == 8
    % Normalize to RPE and RNFL intensities
    % Want to make RPE and RNFL have values of 1
                
    voxel_size = [header.ScaleZ header.ScaleX header.Distance];
    
    if isempty(retina_bds)
        % Estimate RPE and ILM
        [rpe_bdry,ilm_bdry] = quickFindRPE(in_vol,voxel_size);
        rpe_bdry = round(rpe_bdry); 
        ilm_bdry = round(ilm_bdry);
    else
        rpe_bdry = round(retina_bds(:,:,2)); 
        ilm_bdry = round(retina_bds(:,:,1));
    end
                
    % Find the intensity correction factor for the RPE
    minv = 0;
    maxv = round(80/(voxel_size(1)*1000));    
    above = true;
    maxint = true;
    rpe_int_sm = intensityCorrectionFactor(in_vol,rpe_bdry,voxel_size,minv,maxv,above,maxint);
    
    % Protect against very small values
    rpe_int_sm(rpe_int_sm<0.05) = 0.05; 
    
    % Find the intensity correction factor for the RNFL 
    above = false;
    rnfl_int_sm = intensityCorrectionFactor(in_vol,ilm_bdry,voxel_size,minv,maxv,above,maxint);
                
    % Fit linear correction between RPE and RNFL so they are both
    %   approx. 1
    m = (rnfl_int_sm - rpe_int_sm)./(ilm_bdry-rpe_bdry);
    b = rnfl_int_sm - m.*ilm_bdry;

    lin_corr = bsxfun(@plus,bsxfun(@times,shiftdim(m,-1),(1:size(in_vol,1))'),shiftdim(b,-1));

    % Only linear over the retina area
    for i = 1:size(in_vol,2)
        for j = 1:size(in_vol,3)
            lin_corr(1:ilm_bdry(i,j),i,j) = lin_corr(ilm_bdry(i,j),i,j);
            lin_corr(rpe_bdry(i,j):end,i,j) = lin_corr(rpe_bdry(i,j),i,j);
        end
    end                

    out_vol = in_vol./lin_corr*0.8;

    out_vol(out_vol>1) = 1;
    out_vol(out_vol<0) = 0;
    
elseif method == 9
    % Normalize to smoothly varying bias field based on RPE intensity

    voxel_size = [header.ScaleZ header.ScaleX header.Distance];
    
    if isempty(retina_bds)
        % Estimate RPE and ILM
        [rpe_bdry,~] = quickFindRPE(in_vol,voxel_size);
        rpe_bdry = round(rpe_bdry); 
    else
        rpe_bdry = round(retina_bds(:,:,2)); 
    end

    % Find the intensity correction factor for the RPE
    minv = 0;
    maxv = round(80/(voxel_size(1)*1000));    
    above = true;
    maxint = true;
    [~, rpe_int_sm] = intensityCorrectionFactor(in_vol,rpe_bdry,voxel_size,minv,maxv,above,maxint);
    
    % Protect against very small values
    rpe_int_sm(rpe_int_sm<0.05) = 0.05;    
    
    % Divide by correction factor to make RPE intensities approx. 1
    out_vol = bsxfun(@rdivide,in_vol,shiftdim(rpe_int_sm,-1))*.90;
    out_vol(out_vol>1) = 1;    
    
elseif method == 10
    % Normalize to RPE intensity only on B-scans

    voxel_size = [header.ScaleZ header.ScaleX header.Distance];
    
    if isempty(retina_bds)
        % Estimate RPE and ILM
        [rpe_bdry,~] = quickFindRPE(in_vol,voxel_size);
        rpe_bdry = round(rpe_bdry); 
    else
        rpe_bdry = round(retina_bds(:,:,2)); 
    end

    % Find the intensity correction factor for the RPE
    minv = 0;
    maxv = round(80/(voxel_size(1)*1000));    
    above = true;
    maxint = true;
    [~, ~, bscan_fit] = intensityCorrectionFactor(in_vol,rpe_bdry,voxel_size,minv,maxv,above,maxint);
    
    % Protect against very small values
    bscan_fit(bscan_fit<0.05) = 0.05;    
    
    % Divide by correction factor to make RPE intensities approx. 1
    out_vol = bsxfun(@rdivide,in_vol,shiftdim(bscan_fit,-1));
    out_vol(out_vol>1) = 1;
elseif method == 11
    % Normalize by N3O
    [out_vol,~]=n3oct(in_vol,header);
end

% out_vol(out_vol > 1) = 1;
% out_vol(out_vol < 0) = 0;
out_vol(isnan(in_vol)) = nan;

end

