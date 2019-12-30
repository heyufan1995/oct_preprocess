function vol_org = reconstruct_layer(vol,record_params)
% reconstruct from patch and index
% Args:
% vol: output from the network [patch,img_rows,img_cols ,label]
%      or [patch,img_rows,img_cols]
% index: the patch index in the bscan
% prep_record_params: therecord_params returned by Preprocess.m
% Return:
% vol_org: [ori_img_rows, ori_img_cols, Bscans]
disp('Reconstruct layer masks');
% permute vol into [patch, label, img_rows, img_cols]
vol = permute(vol,[1,4,2,3]);
patchindex = record_params.index;
dilation = record_params.dilation;
tv = record_params.tv;
shifts = record_params.shifts;
% sizes of original image: img_rows, img_cols, Bscans.
sizes = record_params.sizes;
img_cols = record_params.options.img_cols;
img_rows = record_params.options.img_rows;

labels = size(vol,2);
segs = length(patchindex);
% standard input shape to graph cut
vol_org = zeros(sizes(1),sizes(2),sizes(3),labels);
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
    for jj = 1:labels
        for kk = 1:segs
                dil = dilation((ii-1)*segs+kk);
                tv_new = tv - (dil-img_rows);
                bv = tv_new + dil-1;
                img = squeeze(vol((ii-1)*segs+kk,jj,:,:));
                if dil~=img_rows
                    [x,y] = meshgrid(1:img_cols,1:img_rows);
                    [xx,yy] = meshgrid(1:img_cols,linspace(1,img_rows,dil)); 
                    img = interp2(x,y,img,xx,yy);
                end
                vol_org(tv_new:bv,patchindex(kk):patchindex(kk)+img_cols-1,ii,jj)=...
                vol_org(tv_new:bv,patchindex(kk):patchindex(kk)+img_cols-1,ii,jj) ...
                +img.*repmat(overlaps(kk,:),[dil,1]);
        end
    end
end
for jj = 1:labels
    vol_org(:,:,:,jj) = retinaFlatten(vol_org(:,:,:,jj),-shifts,'nearest'); 
end
[~,vol_org] = max(vol_org,[],4);
vol_org = vol_org -1;
end

