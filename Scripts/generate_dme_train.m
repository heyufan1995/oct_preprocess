% Generate Training Data
% Results:
    % scan: image data [timepoint, Bscan, img_rows, img_cols]
    % bds: manual boundary [timepoint, Bscan, bdc, img_cols]
    %     the label comes from matlab, remenber -1 in boundary points
    % lesion: manual lesion mask [timepoint, Bscan, img_rows, img_cols]
addpath(genpath('../OCTMatTool'));
addpath(genpath('dme'));
flat_options;
% replace filenames to your own data path
filenames = {    
             'E:\Research\Dataset\OCT\dme\DME\Subject_01.mat'
             'E:\Research\Dataset\OCT\dme\DME\Subject_02.mat' 
             'E:\Research\Dataset\OCT\dme\DME\Subject_03.mat'
             'E:\Research\Dataset\OCT\dme\DME\Subject_04.mat'
             'E:\Research\Dataset\OCT\dme\DME\Subject_05.mat'
             'E:\Research\Dataset\OCT\dme\DME\Subject_06.mat'
             'E:\Research\Dataset\OCT\dme\DME\Subject_07.mat'
             'E:\Research\Dataset\OCT\dme\DME\Subject_08.mat'
             'E:\Research\Dataset\OCT\dme\DME\Subject_09.mat'
             'E:\Research\Dataset\OCT\dme\DME\Subject_10.mat'
             };
mkdir label;mkdir image;
for i = 1:length(filenames)
    filename = filenames{i};
    options.segfile = filename;
    [data, record_params] = Preprocess(filename,options);
    scan = shiftdim(permute(data.flat_vol,[3,1,2]),-1);
    bds = shiftdim(permute(data.bds,[2,3,1]),-1);
    lesion = shiftdim(permute(data.lesion,[3,1,2]),-1);
    % save data
    % crop to the center 500 pixels as the MICCAI 2017 paper did
    % val_clims = 135:634;
    val_clims = 1:768;
    for j = 1:size(scan,2)
        image = squeeze(scan(1,j,:,val_clims));
        label = {};
        label.bds = squeeze(bds(1,j,:,val_clims));
        label.lesion = squeeze(lesion(1,j,:,val_clims));
        label = jsonencode(label);
        [~,name,~] = fileparts(filename);        
        fid = fopen(fullfile('label',sprintf('%s_%d.txt',name,j)),'wt');
        fprintf(fid, label);
        fclose(fid);        
        imwrite(image,fullfile('image',sprintf('%s_%d.png',name,j)),'PNG');    
        fprintf('.')
    end
end


  