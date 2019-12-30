function bd_pts = reconstruct_boundary(bds,record_params)
% reconstruct from patch and index to boundary position
% bds [patch_number, img_cols*boundaries] 
% bd_pts [img_cols,Bscans,boundaries]
disp('Reconstruct boundary');

patchindex = record_params.index;
dilation = record_params.dilation;
tv = record_params.tv;
shifts = record_params.shifts;
% sizes of original image: img_rows, img_cols, Bscans.
sizes = record_params.sizes;
img_cols = record_params.options.img_cols;
img_rows = record_params.options.img_rows;

bdc = size(bds,2)/img_cols;
segs = length(patchindex);
bd_pts = zeros(sizes(2),sizes(3),bdc);

repeatition = zeros(1,sizes(2));
for ii = 1:segs
    repeatition(patchindex(ii):patchindex(ii)+img_cols-1) = repeatition(patchindex(ii):patchindex(ii)+img_cols-1)+1;
end
patchindex_end = find(diff(repeatition)<0);
%% calculate weights for concatenation
overlaps = ones(segs,img_cols);
for ii = 1:segs
    if ii == 1
        right = patchindex_end(ii);
        left = double(patchindex(ii+1));
        x = left:right;
        overlaps(ii,x) = 1 - 1/(right-left)*(x-left);
    elseif ii == segs
        right = patchindex_end(ii-1);
        left = double(patchindex(ii));  
        x = left:right;
        overlaps(ii,x-left+1) = 1/(right-left)*(x-left);
    else
        right1 = patchindex_end(ii-1);
        left1 = double(patchindex(ii));  
        x = left1:right1;
        overlaps(ii,x-left1+1) = 1/(right1-left1)*(x-left1);    
        right2 = patchindex_end(ii);
        left2 = double(patchindex(ii+1));
        x = left2:right2;
        overlaps(ii,x-left1+1) = 1 - 1/(right2-left2)*(x-left2);
    end
end
%% concatenate
for ii = 1:sizes(3)
    for kk = 1:segs
        bd_patch = tv-(dilation((ii-1)*segs+kk)-img_rows)-1 ...
                   +cumsum(reshape(bds((ii-1)*segs+kk,:) ... 
                   *dilation((ii-1)*segs+kk)/img_rows,img_cols,bdc),2);
        bd_pts(patchindex(kk):patchindex(kk)+img_cols-1,ii,:)=...
        squeeze(bd_pts(patchindex(kk):patchindex(kk)+img_cols-1,ii,:))+bd_patch.*repmat(overlaps(kk,:),bdc,1)';
    end
end
bd_pts = bd_pts + repmat(shifts,[1,1,bdc]);
disp('done!');
end

