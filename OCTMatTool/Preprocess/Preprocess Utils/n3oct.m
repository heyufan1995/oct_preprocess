function [oct_vol_n3,gain_field] = n3oct(oct_vol,header,n3params)
% Run N3 for OCT (N3O) intensity normalization and inhomogeneity
% correction, which applies the N3 inhomogeneity correction algorithm to
% OCT data in macular flat space (MFS), thereby allowing better correction
% of individual layers. The imaged area should cover approximately a 6x6 mm
% square around the macula. N3O has been tested on data from both
% Spectralis and Cirrus scanners.
%
% Note that since the macular flat space transformation was trained using
% macular data scanned in a 6x6mm area around the fovea, performance may
% suffer away from this area, particularly near the ONH.
%
% Inputs:
%   oct_vol - NxMxL array (N pixels per A-scan, M A-scans, L B-scans)
%   header - struct with info about the data as output by octReader
%   n3params - struct containing parameters for running N3O. Missing
%     parameters are taken as defaults. Variables include:
%       numIterations - number of iterations of N3 to run (default: 100)
%       controlPointDistance - distance between B-spline control points in
%           the x and y directions for the spline smoothing step; larger
%           values increase the amount of smoothness (units of micron,
%           default: [80 80])
%       lambda - regularization coefficient for the B-spline
%       smoothing step; larger values increase the penalty on the second
%           order derivatives (default: 10000)
%       numHistogramBins - number of histogram bins used in N3 (default:
%           200)
%       convergenceThreshold - convergence threshold of N3; algorithm will
%           stop if this threshold is reached before numIterations are
%           reached (default: 0.001)
%       fwhm - full width half maximum of the gain field distribution in N3
%           (default: 0.1)
%       filterNoise - assumed noise level used by the Wiener filter in N3
%           (default: 0.01)
%       normalize_field - normalize the gain field to have an average value
%           of 1 (default: true)
%       resize_width - width in pixels to resize every B-scan image in the
%           volume to have (default: 256)
%       resize_height - height in pixels to resize every B-scan image in
%           the volume to have (default: 512)
%       pre_norm_data - intensity normalize the data before running N3O;
%           normalization is done using a robust contrast stretching method
%           used in our prior layer segmentation work (default: false)
%       post_norm_data - intensity normalization after running N3O; uses a
%           histogram peak matching method to rescale the intensities in
%           the vitreous and RPE areas (defualt: true)
%       verbose - output algorithm progress to command window (default:
%           true)
%
% For futher algoritm details, see the following publication:
%   A. Lang, A. Carass, B. Jedynak, S. D. Solomon, P. A. Calabresi, J. L.
%       Prince. �Intensity inhomogeneity correction of macular OCT using N3
%       and retinal flatspace.� Proc. International Symposium on Biomedical
%       Imaging (ISBI'16), 197-200, Prague, Czech Republic, April 2016.
%
% Andrew Lang
% 7/13/2016

if nargin < 2
    error('n3oct requires at least 2 inputs (oct_vol and header)')
end
if nargin < 3
    % Use default params
    n3params = [];  
end

n3params = checkOpts(n3params);

% Simple check for Cirrus data
if isa(oct_vol,'uint8')
    oct_vol = im2double(oct_vol);
end

%% Flat space parameters
% These are fixed
flat_params.regression_file = 'flat_space_regression_params_poly2_reg_xydiff_lambda_1060.mat';
flat_params.type = 3;
flat_params.grid_spacing = 4; % average vertical spacing of pixels (micron)
flat_params.pixels_per_layer = 10; % if type == 4
flat_params.pad_above = 60; % Padding above and below ILM and BrM (micron)
flat_params.pad_below = 60;
flat_params.interp_method = 'cubic';
flat_params.smooth_def_interp = true; % smoother deformation field 
flat_params.exact_inverse = false; % only for smooth_def_interp = true
flat_params.verbose = n3params.verbose;

%% Run algorithm

% --- Initial normalization --- %
if n3params.pre_norm_data
    % Note that we can't do this after the downsampling since the output
    % is corrected based on the normalized data
    if n3params.verbose
        fprintf('Initial intensity normalization...\n')
    end
    oct_vol = normalizeOCTVolume(oct_vol,2,header);
end

% --- Downsample data --- %
full_size = size(oct_vol);
% Don't upsample
if n3params.resize_height > full_size(1)
    n3params.resize_height = full_size(1);
end
if n3params.resize_width > full_size(2)
    n3params.resize_width = full_size(2);
end
new_size = [n3params.resize_height n3params.resize_width full_size(3)];
oct_vol_fs = oct_vol; % keep full size vol for output
if ~isequal(full_size,new_size)
    header.ScaleX = header.ScaleX*full_size(2)/new_size(2);
    header.ScaleZ = header.ScaleZ*full_size(1)/new_size(1);
    header.SizeX = new_size(2);
    header.SizeZ = new_size(1);
    oct_vol = imresize(oct_vol,new_size(1:2));
end

% --- Get retina boundaries --- %
if n3params.verbose
    fprintf('Finding retina boundaries...\n')
end
if isfield(header,'scanner_type')
    scanner_type = 'cirrus';
else
    % Simple guess based on number of pixels per A-scan
    if size(oct_vol,1) == 1024
        scanner_type = 'cirrus';
    else
        scanner_type = 'spectralis';
    end
end
[~, ~, bds] = retinaDetector(oct_vol,header,scanner_type);

% --- Convert to flat space --- %
if n3params.verbose
    fprintf('Converting to flat space...\n')
end
flat_params.retina_boundaries = bds;
oct_vol_flat = convertVolToFlatSpace(oct_vol,header,flat_params);

% --- Build template intensity volume --- %
norm_int = sum(sum(oct_vol_flat,2),3)/(size(oct_vol_flat,2)*size(oct_vol_flat,3));
norm_vol_flat = repmat(norm_int,[1 size(oct_vol_flat,2) size(oct_vol_flat,3)]);
norm_vol_flat = norm_vol_flat/max(norm_vol_flat(:));
    
% --- N3 inhomogeneity correction --- %
if n3params.verbose
    fprintf('Estimating gain field...')
end
mask = true(size(oct_vol_flat));
t = tic;
gain_field = zeros(size(oct_vol_flat));
spMat = [];
for ii = 1:size(oct_vol_flat,3)
    if n3params.verbose
        fprintf('B-scan %d of %d\n',ii,size(oct_vol_flat,3));
    end
    [gain_field(:,:,ii),~,spMat] = myn3(oct_vol_flat(:,:,ii),...
                norm_vol_flat(:,:,ii),mask(:,:,ii),header,n3params,spMat);
end
tm = toc(t);
if n3params.verbose
    fprintf('took %f seconds\n',tm)
end

% --- Undo flat space conversion --- %
if n3params.verbose
    fprintf('Unflattening gain field...\n')
end
gain_field = convertVolFromFlatSpace(gain_field,header,flat_params);

% --- Extrapolate gain field --- %
% Smoothly go to 1 outside the flatspace boundaries
if n3params.verbose
    fprintf('Extrapolating gain field...\n')
end
% Extrapolation distance
buf = 150; % micron
buf = round(buf/(header.ScaleZ*1000)); % convert to pixels

gain_field = extrapolateGainField(gain_field,buf);

% --- Upsample result --- %
if ~isequal(full_size,new_size)
    gain_field = imresize(gain_field,full_size(1:2));
end

% --- Corrected volume --- %
oct_vol_n3 = oct_vol_fs./gain_field;

% --- Post-N3O intensity normalization --- %
if n3params.post_norm_data
    if ~isequal(full_size,new_size)
        bds = imresize(bds,full_size(2:3));
        bds = (bds-0.5)*full_size(1)/new_size(1)+0.5;
    end
    % Post-normalize using histogram peak matching
     oct_vol_n3 = normPeaks(oct_vol_n3,header,bds);
end

function gain_field_extrap = extrapolateGainField(gain_field,buf)
% Since the gain field was computed on a masked region, extend the gain
% field away from the top and bottom crop boundaries by linearly
% extrapolating such that the gain field has a value of one at a distance
% of 'buf' pixels and beyond.

sz = size(gain_field);

% Masked areas have a nan value
nv = isnan(gain_field);
    
% Get top and bottom boundary indices
[~,top_inds] = min(nv,[],1);
[~,bot_inds] = min(flipdim(nv,1),[],1);
bot_inds = sz(1)-bot_inds+1;

% Get gain values at the top and bottom boundaries by converting column
% indices to volume indices
[x,y] = meshgrid(1:size(top_inds,3),1:size(top_inds,2));
x = shiftdim(x,-1);
y = shiftdim(y,-1);
top_vals = gain_field(top_inds + (y-1)*sz(1) + (x-1)*sz(1)*sz(2));
bot_vals = gain_field(bot_inds + (y-1)*sz(1) + (x-1)*sz(1)*sz(2));

% Extrapolate these values linearly to 1 over a fixed distance
inds = (1:buf)';
% Set up extrapolation for each A-scan by using values normalized so they
% are between zero and one - this will be scaled and translated to get the
% final extrapolation values
extrap_n = repmat(inds/(buf+1),[1 sz(2),sz(3)]);
% Extrapolated values
extrap_vals_top = bsxfun(@times,extrap_n,(top_vals-1))+1;
extrap_vals_bot = bsxfun(@times,1-extrap_n,(bot_vals-1))+1;
% Column indices of extrapoled data
extrap_inds_top_col = bsxfun(@plus,inds-(buf+1),top_inds);
extrap_inds_bot_col = bsxfun(@plus,inds,bot_inds);

% Convert column to volume indices
extrap_inds_top = bsxfun(@plus,extrap_inds_top_col,...
                               (y-1)*sz(1)+(x-1)*sz(1)*sz(2));
extrap_inds_bot = bsxfun(@plus,extrap_inds_bot_col,...
                               (y-1)*sz(1)+(x-1)*sz(1)*sz(2));
                           
% Remove out of bounds points
extrap_inds_top = extrap_inds_top(extrap_inds_top_col>0);
extrap_vals_top = extrap_vals_top(extrap_inds_top_col>0);
extrap_inds_bot = extrap_inds_bot(extrap_inds_bot_col<=size(gain_field,1));
extrap_vals_bot = extrap_vals_bot(extrap_inds_bot_col<=size(gain_field,1));

% Set extrapolated values in volume
gain_field_extrap = gain_field;    
gain_field_extrap(extrap_inds_top) = extrap_vals_top;
gain_field_extrap(extrap_inds_bot) = extrap_vals_bot;
gain_field_extrap(isnan(gain_field_extrap)) = 1;


function optsOut = checkOpts(opts)
% Function to check user input options with default options. Default
% options are added to missing or empty user options.

defaultOpts = struct('numIterations',50,...
                     'controlPointDistance',[80 80],...
                     'lambda',10000,...
                     'numHistogramBins',200,...
                     'convergenceThresh',0.001,...
                     'fwhm',0.1,...
                     'filterNoise',0.01,...
                     'normalize_field',true,...
                     'resize_width',256,...
                     'resize_height',512,...
                     'pre_norm_data',false,...
                     'post_norm_data',true,...
                     'verbose',true...
                     );

optsOut = defaultOpts;
if ~isempty(opts)
    default_fields = fieldnames(defaultOpts);
    for i=1:length(default_fields)
        if isfield(opts,default_fields{i})
            if ~isempty(opts.(default_fields{i}))
                optsOut.(default_fields{i}) = opts.(default_fields{i});
            end
        end
    end
end