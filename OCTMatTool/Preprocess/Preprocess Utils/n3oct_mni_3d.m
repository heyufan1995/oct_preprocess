function [oct_vol_n3,gain_field] = n3oct_mni_3d(oct_vol,header,n3params)
% Run N3 code on OCT data without using flat space

if nargin < 2
    error('n3oct_mni_3d requires at least 2 inputs (oct_vol and header)')
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

%% Run algorithm

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
if n3params.use_mask || n3params.crop_data
    if n3params.verbose
        fprintf('Finding retina boundaries...\n')
    end    
    retina_mask = retinaDetector(oct_vol,header,scanner_type);
    retina_mask = retina_mask > 0;
else
    retina_mask = [];
end

% --- Initial normalization --- %
if n3params.pre_norm_data
    if n3params.verbose
        fprintf('Initial intensity normalization...\n')
    end
    oct_vol = normalizeOCTVolume(oct_vol,2,header);
end

% --- Crop data to retina limits (with some padding) --- %
if n3params.crop_data    
    pad_size = round(60/(header.ScaleZ*1000));
    tl = find(sum(retina_mask(:,:),2),1,'first') - pad_size; % top limit
    bl = find(sum(retina_mask(:,:),2),1,'last') + pad_size; % bottom limit
    if tl < 1, tl = 1; end
    if bl > size(oct_vol,1), bl = size(oct_vol,1); end

    oct_vol_crop = oct_vol(tl:bl,:,:);
    
    if n3params.use_mask
        retina_mask = retina_mask(tl:bl,:,:);
    else
        retina_mask = [];
    end
else
    oct_vol_crop = oct_vol;
end

% --- N3 inhomogeneity correction --- %
if n3params.verbose
    fprintf('Estimating gain field...')
end
t = tic;
if n3params.use_mask
    mask = retina_mask;
else
    mask = [];
end
gain_field_crop = myn3_3d(oct_vol_crop,oct_vol_crop,mask,header,n3params);
if ~isempty(mask)
    gain_field_crop(~mask) = nan;
end  
tm = toc(t);
if n3params.verbose
    fprintf('took %f seconds\n',tm)
end

% --- Uncrop --- %
if n3params.crop_data
    gain_field = nan(size(oct_vol));
    gain_field(tl:bl,:,:) = gain_field_crop;
else
    gain_field = gain_field_crop;
end

% --- Extrapolate gain field if cropped or masked --- %
if n3params.crop_data || n3params.use_mask
    if n3params.verbose
        fprintf('Extrapolating gain field...\n')
    end
    
    % Extrapolation distance
    buf = 150; % micron
    buf = round(buf/(header.ScaleZ*1000)); % convert to pixels
    
    gain_field = extrapolateGainField(gain_field,buf);
end

% --- Upsample result --- %
if ~isequal(full_size,new_size)
    gain_field = imresize(gain_field,full_size(1:2));
end

% --- Output corrected volume --- %
oct_vol_n3 = oct_vol_fs./gain_field;

% --- Post-N3O intensity normalization --- %
if n3params.post_norm_data
    [~,~,bds] = retinaDetector(oct_vol_n3,header,scanner_type);
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
extrap_inds_top = bsxfun(@plus,(extrap_n-1)*(buf+1),top_inds);
extrap_inds_bot = bsxfun(@plus,extrap_n*(buf+1),bot_inds);

% Convert column to volume indices
extrap_inds_top = bsxfun(@plus,extrap_inds_top,...
                               (y-1)*sz(1)+(x-1)*sz(1)*sz(2));
extrap_inds_bot = bsxfun(@plus,extrap_inds_bot,...
                               (y-1)*sz(1)+(x-1)*sz(1)*sz(2));

% Set extrapolated values in volume
gain_field_extrap = gain_field;    
gain_field_extrap(extrap_inds_top) = extrap_vals_top;
gain_field_extrap(extrap_inds_bot) = extrap_vals_bot;
gain_field_extrap(isnan(gain_field_extrap)) = 1;

function optsOut = checkOpts(opts)
% Function to check user input options with default options. Default
% options are added to missing or empty user options.

% Note controlPointDistance is larger here, run out of memory with smaller
% distances
defaultOpts = struct('numIterations',50,...
                     'controlPointDistance',[80 500 500],...
                     'lambda',10000,...
                     'numHistogramBins',200,...
                     'convergenceThresh',0.001,...
                     'fwhm',0.1,...
                     'filterNoise',0.01,...
                     'normalize_field',false,...
                     'resize_width',256,...
                     'resize_height',512,...
                     'pre_norm_data',false,...
                     'post_norm_data',true,...
                     'verbose',true,...
                     'crop_data',true,...
                     'use_mask',true...
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