function [img_vol,retina_mask,bds,shifts] = preprocessData_Stargart(img_vol,header,pp_params,scanner_type,fids,seg_pts)

% if pp_params.normalize > 0
%     mfprintf(fids,'Normalizing intensities...');
% %     try
%         if pp_params.normalize == 1
%             % Median filter contrast stretching normalization
%             img_vol = normalizeOCTVolume(img_vol,2,header);
%         elseif pp_params.normalize == 2
%             % Attenuation correction
%             img_vol = normalizeOCTVolume(img_vol,5,header);
%         else
%             img_vol = normalizeOCTVolume(img_vol,pp_params.normalize,header);
%         end
%         mfprintf(fids,'done!\n');
% %     catch err
% %         mfprintf(fids,' error!\nError message: %s\n\n',err.message);
% %         img_vol = [];
% %         return        
% %     end
% end
if nargin < 6
    seg_pts = [];
end

%-- Generate retina mask
mfprintf(fids,'Detecting retina boundaries...');
try
    if pp_params.fast_rpe
        if strcmp(scanner_type,'cirrus')
            img_vol = im2double(img_vol);
        end
        [~, ~, shifts, retina_mask, bds, nbpt] = quickFindRPE(img_vol,header);
    else
        if strcmp(scanner_type,'spectralis')
            [retina_mask, shifts, bds, nbpt] = retinaDetector_Stargart(img_vol,header,pp_params.retinadetector_type,false); 
        else
            % Need to median filter first
            sz = size(img_vol);
            if length(sz) == 2
                sz(3) = 1;
            end
            dn_k = [3 3 1];
            img_vol_mf = permute(img_vol,[2 1 3]);
            img_vol_mf = medfilt2(img_vol_mf(:,:),[dn_k(2) dn_k(1)],'symmetric');
            img_vol_mf = reshape(img_vol_mf,sz(2),sz(1),sz(3));
            img_vol_mf = permute(img_vol_mf,[2 1 3]);

            img_vol_mf = im2double(img_vol_mf);
            img_vol = im2double(img_vol);

            [retina_mask, shifts, bds, nbpt] = retinaDetector(img_vol_mf,header,pp_params.retinadetector_type,false); 
        end
    end
    mfprintf(fids,'done! (%d outlier points)\n',nbpt);

    retina_mask = retina_mask > 0;

    if nbpt > 0.5*size(img_vol,2)
        mfprintf(fids,['Warning: poor fit of retina boundaries '...
                       'detected (%d outlier points). Check for '...
                       'artifacts in the data.\n'],nbpt);
    end
catch err
    mfprintf(fids,' error!\nError message: %s\n\n',err.message);
    img_vol = [];
    return
end
if ~isempty(seg_pts)
    % Use manual segmentation
    rpe = seg_pts(:,:,1);
    isos = seg_pts(:,:,7);
    bm = seg_pts(:,:,9);
    
    % Don't overwrite bds since the normalization should use the flattened
    % boundaries
    bds_seg = cat(3,rpe,isos,bm);
    bds_r = round(bds_seg);
%     sz = size(rpe);
    
    if strcmp(scanner_type,'cirrus')
        img_vol = im2double(img_vol);
    end
    retina_mask = false(size(img_vol));
    for i = 1:size(img_vol,2)
        for j = 1:size(img_vol,3)
            if bds_r(i,j,1) > 0
                retina_mask(bds_r(i,j,1):bds_r(i,j,3),i,j) = true;
            end
        end
    end
    
    shifts = bsxfun(@minus,bm,mean(bm,1)+(round(size(img_vol,1)/2)-mean(bm,1)));
end

if pp_params.normalize > 0
    mfprintf(fids,'Normalizing intensities...');
%     try
        if pp_params.normalize == 1
            % Median filter contrast stretching normalization
            img_vol = normalizeOCTVolume(img_vol,2,header);
        elseif pp_params.normalize == 2
            % Attenuation correction
            img_vol = normalizeOCTVolume(img_vol,4,header);
        else
            bds_n = bds(:,:,[1 3]);
            img_vol = normalizeOCTVolume(img_vol,pp_params.normalize,header,bds_n);
        end
        mfprintf(fids,'done!\n');
%     catch err
%         mfprintf(fids,' error!\nError message: %s\n\n',err.message);
%         img_vol = [];
%         return        
%     end
end

%-- Flatten to bottom boundary
if pp_params.flatten == true
    if isfield(pp_params,'flatten_to_isos') && pp_params.flatten_to_isos
        if isempty(seg_pts)
            isos = bds(:,:,2);
        else
            isos = bds_seg(:,:,2);
        end
        shifts = bsxfun(@minus,isos,mean(isos,1)+(round(size(img_vol,1)/2)-mean(isos,1)));
    end
    
    % Make sure the data is still within the image after flattening
    tb = bds(:,:,1)-shifts;
    if any(tb < 0)
        shifts = shifts + min(tb(:));
        % Center
        shifts = shifts - min(min((size(img_vol,1)-(bds(:,:,end)-shifts))))/2;
    end
    
    mfprintf(fids,'Flattening data...');
    try
        img_vol = retinaFlatten(img_vol,shifts,'linear');
        retina_mask = retinaFlatten(retina_mask,shifts,'nearest');
        mfprintf(fids,'done!\n');
    catch err
        mfprintf(fids,' error!\nError message: %s\n\n',err.message);
        img_vol = [];
        return
    end
end