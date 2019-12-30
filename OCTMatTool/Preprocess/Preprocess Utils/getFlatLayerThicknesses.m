function N = getFlatLayerThicknesses(flat_params,reg_boundaries)
% Given the regression boundaries, determine how many pixels to have in
% flat space for each layer
%
% The number of pixels depends on the values of flat_params.type and
% flat_params.grid_spacing. The variable flat_params.type can take values:
%   1 - Layer thicknesses equal to average thickness of each layer in
%       reg_boundaries
%   2 - Layer thicknesses equal to the average thickness of each layer from
%       a normative set of data (found in the file flat_space_norm_th.mat,
%       see training/saveFlatSpaceNormativeThickness.m)
%   3 - Same as 2 except thicknesses are equal to the average thickness
%       plus one standard deviation (thus it each layer is slightly
%       thicker to work better for thicker retinas)
%   4 - Each layer gets the number of pixels given by the parameter
%       flat_params.pixels_per_layer, which can be a scalar applied to each
%       layer, or a vector with the size of each layer
%
% For types 1, 2, and 3, the size is divided by flat_params.grid_spacing to
% get the number of pixels in each layer (reg_boundaries and the values in
% flat_space_norm_th.mat are in units of micron, so
% flat_params.grid_spacing should have units of pixels/micron)

% Fill any missing options with defaults
default_params = struct('pixels_per_layer',10,'type',2,'grid_spacing',4,...
                        'pad_above',60,'pad_below',60);
flat_params = checkOpts(default_params,flat_params);

N = flat_params.pixels_per_layer;
flat_type = flat_params.type;
grid_spacing = flat_params.grid_spacing;
dist_above = flat_params.pad_above;
dist_below = flat_params.pad_below;

% Regression thicknesses
reg_th = diff(reg_boundaries,1,3);
if flat_type ~= 4
    if flat_type == 1
        % Get number of nodes from average thickness of each layer
        Nw = zeros(size(reg_th,3),1);
        for j = 1:size(reg_th,3)
            % Get average thickness of layer
            lw = reg_th(:,:,j);
            Nw(j) = mean(lw(:));
        end
    else
        % Load average thickness of each layer from file 
        S = load('flat_space_norm_th.mat');
        if flat_type == 2                       
            Nw = S.th_stats.mean;
        else % flat_type == 3
            Nw = S.th_stats.mean + S.th_stats.std;
        end
    end
    % Number of nodes per layer
    N = round(Nw/grid_spacing);
    
    N_above = round(dist_above/grid_spacing);
    N_below = round(dist_below/grid_spacing);
elseif flat_type == 4
    if ~(length(N) == 1 || length(N) == size(reg_th,3) || length(N) == (size(reg_th,3)+2))
        error('Invalid input for number of points per layer in flat space')
    end
    
    if length(N) == 1
        % Input is a scalar for number of points in each layer
        N = repmat(N,[size(reg_th,3) 1]);
    else
        % Input has the number of points for each individual layer
        
    end
    
    N_above = N(1);
    N_below = N(end);
            
    if length(N) == (size(reg_th,3) + 2)
        % Input gives the number of points for each layer plus above and
        % below the retina
        N = N(2:end-1);
    end
else
    error('Invalid flat space type! flat_params.type can take values of 1, 2, or 3.')
end

N = cat(1,N_above,N,N_below);

%%
function optsOut = checkOpts(defaultOpts,opts)
% Function to check user input options with default options. Default
% options are added to missing or empty user options.

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