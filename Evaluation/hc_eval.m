% The surface position groundtruth and prediction from the network 
% starts from 0 (python indexing), but the shortest path results calculated in
% matlab start from 1 and minus 1 is needed. When converting boundary into
% layer masks, plus 1 is needed.

addpath(genpath('../OCTMatTool'))
% path to the result
paths = './hc/';
% original test data name
test_file ={
'hc01_spectralis_macula_v1_s1_R'
'hc02_spectralis_macula_v1_s1_R'
'hc03_spectralis_macula_v1_s1_R'
'hc04_spectralis_macula_v1_s1_R'
'hc05_spectralis_macula_v1_s1_R'
'hc06_spectralis_macula_v1_s1_R'
'hc07_spectralis_macula_v1_s1_R'
'hc08_spectralis_macula_v1_s1_R'
'ms01_spectralis_macula_v1_s1_R'
'ms02_spectralis_macula_v1_s1_R'
'ms03_spectralis_macula_v1_s1_R'
'ms04_spectralis_macula_v1_s1_R'
'ms05_spectralis_macula_v1_s1_R'
'ms06_spectralis_macula_v1_s1_R'
'ms07_spectralis_macula_v1_s1_R'
'ms08_spectralis_macula_v1_s1_R'
'ms09_spectralis_macula_v1_s1_R'
'ms10_spectralis_macula_v1_s1_R'
'ms11_spectralis_macula_v1_s1_R'
'ms12_spectralis_macula_v1_s1_R'   
    };
% visualize
dplot = false;
% evaluation metric
eval_dice = false;
% evaluate shortest path
sp = false;
if sp
    fprintf('calculating shortest path\n')
end
% dice from the mask branch
dice_dl = zeros(8,length(test_file));
% dice from the boundary branch
dice_dl_r = zeros(8,length(test_file));
% mean absolute distance for each bscan
diff_dl_b = zeros(9,49,length(test_file));
% mean absolute distance for each volume
diff_dl = zeros(9,length(test_file));
% mean signed distance for each volume
msd_dl = zeros(9,length(test_file));
% rooted mean square distance for each volume
rm_dl = zeros(9,length(test_file));

for i = 1:length(test_file)
    fprintf([num2str(i),'\n']);
    filename = test_file{i};
    % H, W, B
    mask_pred = zeros(128, 1024, 49);
    mask_gt = zeros(128, 1024, 49);
    % W, B ,C
    bds_pred = zeros(1024, 49, 9);
    bds_gt = zeros(1024, 49, 9);
    for idx = 1:49
        bscan = [filename, '_', num2str(idx), '_mean'];
        matpath = [paths,bscan];
        s = load(matpath);
        mask_gt(:,:,idx) = s.mask_gt;
        [~, mp] = max(s.mask_pred); 
        mask_pred(:,:,idx) = squeeze(mp) - 1;        
        if sp
            if ~isfield(s,'bds_pred_sp')
                tic
                s.bds_pred_sp = dijkstra_oct(s);
                toc
                save(matpath,'-struct', 's');
            end            
            bds_pred(:,idx,:) = permute(s.bds_pred_sp,[3,1,2]) - 1;           
        else       
            bds_pred(:,idx,:) = permute(s.bds_pred,[3,1,2]);
        end        
        bds_gt(:,idx,:) = permute(s.bds_gt,[2,3,1]);
        if dplot
            figure(1);
            imagesc(s.img);colormap gray;axis image;axis off;
            line(1:1024,squeeze(s.bds_pred)+1,'LineWidth',2);
            figure(2);
            imagesc(s.img);colormap gray;axis image;axis off;
            line(1:1024,squeeze(s.bds_pred_sp),'LineWidth',2);            
            figure(3);
            imagesc(s.img);colormap gray;axis image;axis off;
            line(1:1024,squeeze(s.bds_gt)+1,'LineWidth',2);            
            pause;
        end
            
    end
    % evaluation
    % dice
    if eval_dice
        for ii = 1:8
            l1 = mask_gt(:)==ii;
            l2 = mask_pred(:)==ii;
            dice_dl(ii,i) = 2*sum(l1.*l2)/(sum(l1)+sum(l2)); 
            layerr = convertBoundariesToLabels(bds_pred+1,[128,1024,49]);
            l2 = layerr(:)==ii;
            dice_dl_r(ii,i) = 2*sum(l1.*l2)/(sum(l1)+sum(l2));           
        end
    end
    % evaluation 
    % mae
    diff_dl_b(:,:,i) = squeeze(mean(abs(bds_gt-bds_pred),1))'*3.9;
    diff_dl(:,i) = mean(mean(abs(bds_gt-bds_pred),1),2)*3.9;
    % msd
    msd_dl(:,i) = mean(mean((bds_gt-bds_pred),1),2)*3.9;
    % rmse
    rm_dl(:,i) = squeeze(sqrt(mean(mean((bds_gt-bds_pred).^2,2),1)))*3.9;
end