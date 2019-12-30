function boundary_points = convertBoundariesFromFlatSpace_new2(boundary_points_flat,header,flat_params,reg_boundaries)
% method = 1 - flatten to average boundary positions
% method = 2 - flatten to regression learned boundary positions
%
% To add? Flatten to input boundary positions (flatten to each boundary
% given and linearly interpolate between boundaries)

if nargin < 4
    reg_boundaries = [];
end

%% Parameters

grid_spacing = flat_params.grid_spacing;

sm = flat_params.smooth_def_interp;

%% Get estimated position of each boundary from regression
% These are the flat space boundaries

if isempty(reg_boundaries)
    % Compute regression boundaries
    reg_boundaries = thicknessRegression(flat_params,header);

    % Thickness of regression layers
    N = getFlatLayerThicknesses(flat_params,reg_boundaries);
else
    N = flat_params.pixels_per_layer;
end
N_above = N(1);
N_below = N(end);

% Number of points to add for each layer
N_csum = cumsum(cat(1,0,N(1:end-1))) + (0:(length(N)-1))';

%% Undo flat space

% Get flat space boundary indices/positions
fs_boundaries = N_csum(2:end);
fs_boundaries = [1; fs_boundaries; fs_boundaries(end)+N_below];

% Add a boundary above and below the retina which represent the top and
% bottom of flat space
reg_boundaries = cat(3,reg_boundaries(:,:,1)-grid_spacing*N_above,reg_boundaries);
reg_boundaries = cat(3,reg_boundaries,reg_boundaries(:,:,end)+grid_spacing*N_below);

% Interpolate at given boundary points
reg_boundaries = permute(reg_boundaries,[3 1 2]);
boundary_points_flat = permute(boundary_points_flat,[3 1 2]);
if ~sm
    % Start by getting the linear coefficients
    %   note: each linear piece takes the form y = a*(x-x_c) + b, where x_c
    %   is the left control point of that piece
    pp = interp1(fs_boundaries,reg_boundaries(:,:),'linear','pp');
    pp.coefs = reshape(pp.coefs,[size(reg_boundaries,2),size(reg_boundaries,3),length(fs_boundaries)-1 2]);
    
    % Linear extrapolation - since the parameterization is centered on the
    % left control point, the extrapolated slopes are the same as the end
    % pieces, but the intercept is shifted
    pp.breaks = [pp.breaks(1)-1 pp.breaks pp.breaks(end)+1];
    % Slope (doesn't change)
    m1 = pp.coefs(:,:,1,1);
    m2 = pp.coefs(:,:,end,1);    
    % Intercept
    b1 = pp.coefs(:,:,1,2) - m1.*(pp.breaks(2)-pp.breaks(1));
    b2 = pp.coefs(:,:,end,2) - m2.*(pp.breaks(end-2)-pp.breaks(end-1));    
    % Concatenate to pp structure
    ppc1 = permute(cat(3,m1,b1),[1 2 4 3]);
    ppc2 = permute(cat(3,m2,b2),[1 2 4 3]);
    pp.coefs = cat(3,ppc1,pp.coefs,ppc2);
    pp.pieces = pp.pieces+2;
    
    % Bin the data to get the section the data lies in
    indices = ones(size(boundary_points_flat),'uint8');
    for ii = 1:length(fs_boundaries)
        c = boundary_points_flat > fs_boundaries(ii);
        indices(c) = ii+1;
    end

    % Polynomial equation for corresponding piece
    c_inds = zeros([size(boundary_points_flat),2]);
    for i = 1:size(boundary_points_flat,2)
        for j = 1:size(boundary_points_flat,3)
            c_inds(:,i,j,:) = pp.coefs(i,j,indices(:,i,j),:);
        end
    end
    
    % Evaluate polynomial
    boundary_points = c_inds(:,:,:,1).*(boundary_points_flat - pp.breaks(indices)) + ...
                      c_inds(:,:,:,2);
    
    
%     for i = 1:size(boundary_points_flat,2)
%         for j = 1:size(boundary_points_flat,3)
%             F = griddedInterpolant(fs_boundaries,reg_boundaries(:,i,j),'linear');
% 
%             F.ExtrapolationMethod = 'linear';
%             boundary_points(:,i,j) = F(boundary_points_flat(:,i,j));
%         end
%     end
else
    % Start by getting the cubic polynomial coefficients
    pp = pchip_fs(fs_boundaries',reg_boundaries(:,:)');
    pp.coefs = reshape(pp.coefs,[size(reg_boundaries,2),size(reg_boundaries,3),length(fs_boundaries)-1 4]);
    
    % Linear extrapolation by fitting a line to the end points
    pp.breaks = [pp.breaks(1)-1 pp.breaks pp.breaks(end)+1];
    % Slope
    m1 = (reg_boundaries(2,:,:)-reg_boundaries(1,:,:))/(fs_boundaries(2)-fs_boundaries(1));
    m2 = (reg_boundaries(end,:,:)-reg_boundaries(end-1,:,:))/(fs_boundaries(end)-fs_boundaries(end-1));
    % Intercept
    b1 = reg_boundaries(1,:,:)-m1.*fs_boundaries(1)+m1.*pp.breaks(1);
    b2 = reg_boundaries(end,:,:)-m2.*fs_boundaries(end)+m2.*pp.breaks(end-1);
    
    % Concatenate to pp structure
    ppc1 = cat(1,zeros(size(m1)),zeros(size(m1)),m1,b1);
    ppc2 = cat(1,zeros(size(m2)),zeros(size(m2)),m2,b2);
    ppc1 = permute(ppc1,[2 3 4 1]);
    ppc2 = permute(ppc2,[2 3 4 1]);    
    pp.coefs = cat(3,ppc1,pp.coefs,ppc2);
    pp.pieces = pp.pieces+2;
    
    % Bin the data to get the section the data lies in
    indices = ones(size(boundary_points_flat),'uint8');
    for ii = 1:length(fs_boundaries)
        c = boundary_points_flat > fs_boundaries(ii);
        indices(c) = ii+1;
    end

    % Polynomial equation for corresponding piece
    c_inds = zeros([size(boundary_points_flat),4]);
    for i = 1:size(boundary_points_flat,2)
        for j = 1:size(boundary_points_flat,3)
            c_inds(:,i,j,:) = pp.coefs(i,j,indices(:,i,j),:);
        end
    end
    
    % Polynomial is evaluted relative to the break point position for the
    % section the point is in so we offset the points
    boundary_points_flat = boundary_points_flat - pp.breaks(indices);
    % Free up memory, cirrus data can use a lot
    clear c pp indices    
    % Evaluate polynomial (note: x.*x.*x is much faster than x.^3)
    boundary_points = c_inds(:,:,:,1).*boundary_points_flat.*...
                                       boundary_points_flat.*...
                                       boundary_points_flat + ...
                      c_inds(:,:,:,2).*boundary_points_flat.^2 + ...
                      c_inds(:,:,:,3).*boundary_points_flat + ...
                      c_inds(:,:,:,4);
end

boundary_points = permute(boundary_points,[2 3 1]);

boundary_points = boundary_points/(header.ScaleZ*1000);

