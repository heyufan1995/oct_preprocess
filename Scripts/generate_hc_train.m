function generate_hc_train(filelist,seglist,result)
% Generate Training Data
% train_vol
% Inputs:
    % filelist: txt file with scan file path
    % seglist: txt file with manual segmentation path
    % result: path to save results
% Results:
    % scan: image data [timepoint, Bscan, img_rows, img_cols]
    % bds: manual boundary [timepoint, Bscan, bdc, img_cols]
    %     the label comes from matlab, remenber -1 in boundary points
    % lesion: manual lesion mask [timepoint, Bscan, img_rows, img_cols]
    
% generate_hc_train('./hc/filename.txt','./hc/segname.txt')
addpath(genpath('../OCTMatTool'));
addpath(genpath('hc'));
if nargin < 3
    result = '.';
end
flat_options;
[~, ~, ext] = fileparts(filelist);
if strcmp(ext,'.txt')
    fid = fopen(filelist,'r');
    filenames_list = textscan(fid,'%s','Delimiter','\n');
    filenames_list = filenames_list{1};
    fclose(fid);
end
[~, ~, ext] = fileparts(seglist);
if strcmp(ext,'.txt')
    fid = fopen(seglist,'r');
    segnames_list = textscan(fid,'%s','Delimiter','\n');
    segnames_list = segnames_list{1};
    fclose(fid);
end
mkdir(fullfile(result,'label'));
mkdir(fullfile(result,'image'));
for i = 1:length(filenames_list)
    options.segfile = segnames_list{i};
    filename = filenames_list{i};
    [data, record_params] = Preprocess(filename,options);
    scan = shiftdim(permute(data.flat_vol,[3,1,2]),-1);
    bds = shiftdim(permute(data.bds,[2,3,1]),-1);
    % lesion = shiftdim(permute(data.lesion,[3,1,2]),-1);
    % save data
    for j = 1:size(scan,2)
        image = squeeze(scan(1,j,:,:));
        label = {};
        label.bds = squeeze(bds(1,j,:,:));
        %label.lesion = squeeze(lesion(1,j,:,:));
        label = jsonencode(label);
        [~,name,~] = fileparts(filename);
        fid = fopen(fullfile(result,'label',sprintf('%s_%d.txt',name,j)),'wt');
        fprintf(fid, label);
        fclose(fid);        
        imwrite(image,fullfile(result,'image',sprintf('%s_%d.png',name,j)),'PNG');    
        fprintf('.')
    end
end


  