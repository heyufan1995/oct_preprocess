function reg_boundaries = thicknessRegression(flat_params,header)
% Estimate retina boundaries from the outer retinal boundaries using a
% learned regression defined by the regression file in the flat_params
% input

graph_file = flat_params.regression_file;

% Get position of fovea (to center the graph)
foveaPos = foveaFinder(flat_params.retina_boundaries,header);

retina_boundaries = flat_params.retina_boundaries*header.ScaleZ*1000;

% Load regression parameters
try
    load(graph_file,'reg_params','options');
    reg_params = permute(reg_params,[1 2 4 3]);
catch
    fprintf('Unable to load graph file!\n')
    return
end

% Graph center
i_cen = round([size(reg_params,1) size(reg_params,2)]/2);

if ~strncmp(header.ScanPosition,'OD',2)
    reg_params = flipdim(reg_params,1);    
    i_cen = size(reg_params(:,:,1,1)) - i_cen + 1;
end

% Rotate if necessary
if isfield(header,'angle') && abs(header.angle) > 0.01
    % Rotate and flip to align with fundus
    % Assuming a horizontal scan, data is acquired left to right then
    % bottom to top, so need to permute and flip y axis to go top to bottom
    % then right to left)
    reg_params = flipdim(permute(reg_params,[2 1 3 4]),1);
    scaleX = 6/1024*1000; % training data is 6x6mm, 1024x49 pixels
    scaleY = 6/49*1000;
    i_cen2 = [i_cen(1) size(reg_params,1)-i_cen(2)+1];
    
    [X,Y] = meshgrid(((1:size(reg_params,2))-i_cen2(1))*scaleX,...
                     ((1:size(reg_params,1))-i_cen2(2))*scaleY);
    [Xr,Yr] = rot2d(X,Y,header.angle);
    for i = 1:size(reg_params,3)
        for j = 1:size(reg_params,4)
            reg_params(:,:,i,j) = interp2(X,Y,reg_params(:,:,i,j),Xr,Yr);
        end
    end
    
    % Undo rotate and flip
    reg_params = permute(flipdim(reg_params,1),[2 1 3 4]);
end

% Coordinates of the A-scans in um
sz = size(retina_boundaries);
x = ((1:sz(1))-foveaPos(1))*header.ScaleX;
y = ((1:sz(2))-foveaPos(2))*header.Distance;

% Convert coordinates to graph data space (created from 6x6 mm, 1024x49 px
% data) - rounding is equivalent to nearest neighbor interpolation
xinds = round(x/3*1024/2);
yinds = round(y/3*49/2);

% Align image indices to graph center
xinds = round(xinds + i_cen(1));
yinds = round(yinds + i_cen(2));

% Extract subject/retina aligned parameters
reg_params = reg_params(xinds,yinds,:,:);

% Get node positions for each boundary from retina thickness
ilm = retina_boundaries(:,:,1);
brm = retina_boundaries(:,:,end);
th = brm-ilm;
th_pos = zeros(size(reg_params,1),size(reg_params,2),size(reg_params,3));
for j = 1:size(reg_params,3)
    for k = 0:options.poly_order
        th_pos(:,:,j) = th_pos(:,:,j) + reg_params(:,:,j,k+1).*th.^k;
    end
end

% Get estimated boundary positions in native space

% ILM boundary
ilm = th_pos(:,:,end) + retina_boundaries(:,:,1);
th_pos = cumsum(th_pos(:,:,1:8),3);

reg_boundaries = cat(3,ilm,bsxfun(@plus,ilm,th_pos));

% Make sure points are monotonic by checking derivative along A-scans
nonmono = cat(3,false(size(reg_boundaries(:,:,1))),...
                diff(reg_boundaries,1,3)<=0);
iters = 0;
while nnz(nonmono) > 0
    % Get indices of non-monotonic points
    inds = find(nonmono);
    
    % Set boundary value slightly above the one above it (1 micron above)
    reg_boundaries(inds) = reg_boundaries(inds-size(nonmono,1)*size(nonmono,2))+1;
    
    % Make sure boundaries now monotonic
    nonmono = cat(3,false(size(reg_boundaries(:,:,1))),...
                    diff(reg_boundaries,1,3)<=0);
    iters = iters + 1;
    if iters > 10
        % Avoid endless loop possibility
        warning('Flatspace boundaries are non-monotonic! Could not fix...')
        break;
    end
end
