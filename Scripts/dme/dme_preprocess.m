function [ img_vol ] = dme_preprocess( img_vol )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

img_vol(isnan(img_vol)) = 0;

% Fill in from the left side
[~,inds] = max(img_vol<255,[],2); % Index of first nonzero value
for i = 1:size(img_vol,1)
    for j = 1:size(img_vol,3)
        p = inds(i,1,j);
        if p > 1 && p < i
            if p < (size(img_vol,2)-2)
                % Avoid using low intensity edge pixels
                img_vol(i,1:(p+1),j) = 0; 
            else
                img_vol(i,1:(p-1),j) = 0;
            end
        end
    end
end
% Fill in from the right side
[~,inds] = max(flipdim(img_vol<255,2),[],2); % Index of last nonzero value
inds = size(img_vol,2) - inds + 1;
for i = 1:size(img_vol,1)
    for j = 1:size(img_vol,3)
        p = inds(i,1,j);
        if p < size(img_vol,2) && size(img_vol,2)-p < i
            if p > 2
                % Avoid using low intensity edge pixels
                img_vol(i,(p-1):end,j) = 0; 
            else
                img_vol(i,(p+1):end,j) = 0;
            end
        end
    end
end
% Fill in from the top
mv = mean(img_vol(:));
[~,inds] = max(img_vol<255,[],1); % Index of first nonzero value
for i = 1:size(img_vol,2)
    for j = 1:size(img_vol,3)
        p = inds(1,i,j);
        if p > 1
            if p < (size(img_vol,1)-2)
                % Avoid using low intensity edge pixels
                if img_vol(p+2,i,j) < mv
                    img_vol(1:(p+1),i,j) = 0;
                else
                    % Probably cut through the retina so keep a gradient
                    % here
                    img_vol(1:(p+1),i,j) = 0;
                end
            else
                img_vol(1:(p-1),i,j) = 0;
            end
        end
    end
end
% Fill in from the bottom
[~,inds] = max(flipdim(img_vol<255,1),[],1); % Index of first nonzero value
inds = size(img_vol,1) - inds + 1;
for i = 1:size(img_vol,2)
    for j = 1:size(img_vol,3)
        p = inds(1,i,j);
        if p < size(img_vol,1)
            if p > 2
                % Avoid using low intensity edge pixels
                img_vol((p-1):end,i,j) = 0; 
            else
                img_vol((p+1):end,i,j) = 0;
            end
        end
    end
end
end

