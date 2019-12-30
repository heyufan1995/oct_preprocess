% Default options
% The preprocessing params. use default.
options.preproc_params.normalize = 2;
options.preproc_params.filter = false;
options.preproc_params.filter_kernel = [1 1 3];
options.preproc_params.flatten = 1;
options.preproc_params.fast_rpe = false;
% data options
options.types = 'hc';
options.img_rows = 128;
options.segfile='';
% crop params
options.crop = false;
% other options
options.org = false; % return original img_vol
options.flatvol = true; % return the flattened bscans
options.dplot = true;