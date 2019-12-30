addpath(genpath('../OCTMatTool'))
% path to the result
paths = './dme';
result = dir(fullfile(paths,'*.mat'));
% image sizes
slice = 55;
img_cols = 768;
bds_num = 8;
img_rows = 224;
% combine the results
if isfile(fullfile(paths,'combined.mat'))
    load(fullfile(paths,'combined'))
else
    dlpred = struct();
    dlpred.pred_avg = nan*zeros(slice,img_cols,bds_num);
    dlpred.pred_mask = nan*zeros(slice,img_rows,img_cols);
    dlpred.gt_bds = zeros(slice,img_cols,bds_num);
    dlpred.gt_mask = zeros(slice,img_rows,img_cols);
    dlpred.img = zeros(slice,img_rows,img_cols);
    for i = 1:length(result)
        S=load(fullfile(result(i).folder,result(i).name));
        % image
        dlpred.img(i,:,:) = S.img;
        % ground truth
        dlpred.gt_bds(i,:,:) = squeeze(S.bds_gt)' + 1;
        dlpred.gt_mask(i,:,:) = S.mask_gt;
        % mask prediction
        mask_pred = S.mask_pred;
        [~,mask_pred] = max(double(mask_pred),[],1);
        mask_pred = squeeze(mask_pred);
        dlpred.pred_mask(i,:,:) = mask_pred;
        % coordinates avg
        dlpred.pred_avg(i,:,:) = squeeze(S.bds_pred)' + 1;
    end
    save(fullfile(paths,'combined'),'dlpred');
end
% our pipeline may output an extra surface with no corresponding manual
if bds_num>8
    dlpred.gt_bds(:,:,6) = [];
    dlpred.pred_avg(:,:,6) = [];
end
% baselines
start = 1;
vis = false;
% only evaluate the mid 500 pixels
% (Rathke, et al. "Locally adaptive probabilistic models on pathological oct scans" Miccai 2017)
val_clims = 135:634;
chiu_result = './dme_baseline/chiu';
c_r = load(chiu_result);
chiu_r = c_r.ChiuPred(start:start+54,:,:);
% drop surface values where chiu's results dropped.
gt = dlpred.gt_bds;
temp = nan*ones(size(gt)); 
temp(:,val_clims,:) = gt(:,val_clims,:);
gt = temp;
gt(isnan(chiu_r)) = nan; 
dlpred.pred_avg(isnan(gt)) = nan;
he = squeeze(nanmean(abs(dlpred.pred_avg-gt),2))*3.87;
he_dice = 2*sum((dlpred.pred_mask(:)==10).*(dlpred.gt_mask(:)==9))...
          /(sum(dlpred.pred_mask(:)==10)+sum(dlpred.gt_mask(:)==9));
