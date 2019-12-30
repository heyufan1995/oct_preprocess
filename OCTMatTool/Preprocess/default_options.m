% Default options
% The preprocessing params. use default.
options.preproc_params.normalize = 2;
options.preproc_params.filter = false;
options.preproc_params.filter_kernel = [1 1 3];
options.preproc_params.flatten = 1;
options.preproc_params.fast_rpe = false;
% crop params
options.types = 'hc';
options.img_rows = 128;
options.img_cols = 128;
options.segfile = ''; % label path
options.train = false; % the options.segs can be large if true otherwise not
options.segs = 10; % default patches per bscan
options.crop = true; % always set to true 
% other options
options.resize='';
options.dplot = false;
options.org = false; % return original img_vol
options.flatvol = true; % return the flattened bscans
