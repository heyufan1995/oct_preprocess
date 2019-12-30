% run the script and test if the toolbox is working
% the example flattens the image and crop into patches and then reconstruct
% back to the original image.
addpath(genpath('Reconstruct'));
addpath(genpath('Visualize'));
addpath(genpath('Preprocess'));
default_options;
options.segfile = 'example.mat';
options.org = true;
filename = 'example.vol';
% Test preprocess and reconstruct
[data, record_params] = Preprocess(filename,options);
vol_org = reconstruct_layer(data.Y_train_lesion,record_params);
data.Y_train_b = reshape(data.Y_train_b, size(data.Y_train_b,1),size(data.Y_train_b,2)*size(data.Y_train_b,3));
bd_pts = reconstruct_boundary(data.Y_train_b,record_params);
% Test visualization
vis_opts.pause = true;
vis_opts.save = false;
visdata.img = data.img_vol;
visdata.mask = vol_org;
visdata.bds = bd_pts;
visualizer(visdata, vis_opts);
