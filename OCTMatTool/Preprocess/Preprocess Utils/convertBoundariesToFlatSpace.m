function [boundary_points_flat,flat_boundary_inds] = convertBoundariesToFlatSpace(boundary_points,header,flat_params)
% Convert boundary points in native space to boundary points in flat space.
% Boundary points is an Na x Nb x Np matrix, where Na is the number of
% A-scans, Nb is the number of B-scans, and Np is the number of boundaries.
% Note that Np can also be the number of pixels per A-scan, and thus, can
% also be used to transform the pixel coordinates of the entire volume
% instead of just the boundaries.

% -- Parameters -- %

grid_spacing = flat_params.grid_spacing;
sm = flat_params.smooth_def_interp;

% -- Get position of each boundary in flat space from regression -- %

% Convert to micron
boundary_points = boundary_points*header.ScaleZ*1000;

reg_boundaries = thicknessRegression(flat_params,header);

% -- Get size of flat space -- %

N = getFlatLayerThicknesses(flat_params,reg_boundaries);
N_above = N(1);
N_below = N(end);

% Number of points to add for each layer
N_csum = cumsum(cat(1,0,N(1:end-1))) + (0:(length(N)-1))';

% -- Build flat space -- %

% Get flat space boundary indices/positions
fs_boundaries = N_csum(2:end);
fs_boundaries = [1; fs_boundaries; fs_boundaries(end)+N_below];

% Add a boundary above and below the retina which represent the top and
% bottom of flat space
reg_boundaries = cat(3,reg_boundaries(:,:,1)-grid_spacing*N_above,reg_boundaries);
reg_boundaries = cat(3,reg_boundaries,reg_boundaries(:,:,end)+grid_spacing*N_below);
reg_boundaries = permute(reg_boundaries,[3 1 2]);

% Interpolate deformation field at given points - the tranformation is
% defined from flat space to native space, so we need to first get this
% transformation and then invert it; the inversion is analytical since the
% forward transformation is defined as either piecewise linear
% (smooth_def_interp = false) or cubic hermite spline (smooth_def_interp =
% true)
boundary_points = permute(boundary_points,[3 1 2]);
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
    indices = ones(size(boundary_points),'uint8');
    for ii = 1:length(fs_boundaries)
        c = bsxfun(@minus,boundary_points,reg_boundaries(ii,:,:)) > 0;
        indices(c) = ii+1;
    end
    % Linear equation for corresponding piece
    c_inds = zeros([size(boundary_points),2]);
    for i = 1:size(boundary_points,2)
        for j = 1:size(boundary_points,3)
            c_inds(:,i,j,:) = pp.coefs(i,j,indices(:,i,j),:);
        end
    end
    
    % Solve inverse of each section: y = a*(x-x0) + b for x    
    def_interp = (boundary_points-c_inds(:,:,:,2))./c_inds(:,:,:,1) + pp.breaks(indices);
    
    def_interp = permute(def_interp,[2 3 1]);
else
    % Start by getting the cubic polynomial coefficients
%     pp = pchip(fs_boundaries',reg_boundaries(:,:)');
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
    indices = ones(size(boundary_points),'uint8');
    for ii = 1:length(fs_boundaries)
        c = bsxfun(@minus,boundary_points,reg_boundaries(ii,:,:)) > 0;
        indices(c) = ii+1;
    end
    % Polynomial equation for corresponding piece
    c_inds = zeros([size(boundary_points),4]);
    for i = 1:size(boundary_points,2)
        for j = 1:size(boundary_points,3)
            c_inds(:,i,j,:) = pp.coefs(i,j,indices(:,i,j),:);
        end
    end
    
    % Solve the equation f(x) = y_i for x, which just finding the
    % roots of a cubic polynomial f(x)-y_i
    c_inds(:,:,:,4) = c_inds(:,:,:,4)-boundary_points;
    clear boundary_points c
    c_inds = reshape(c_inds,[],4);
    % Chunk the data to handle large scans (like Cirrus data)
    p = round(linspace(1,size(c_inds,1)+1,20));
    r3 = zeros(size(c_inds,1),3);
    linear_inds = false(size(c_inds,1),1);
    for i = 1:(length(p)-1)    
        [r3p,linear_inds_p] = cs_inv2(c_inds(p(i):(p(i+1)-1),:)); % get the roots
        r3(p(i):(p(i+1)-1),:) = r3p;
        linear_inds(p(i):(p(i+1)-1)) = linear_inds_p;
    end
    clear c_inds
    
    % Extract the real roots that lie in the interval of interest
    % (closest positive root to 0 since each cubic equation is
    % offset by the breakpoints)
    r3lin = r3(linear_inds,1);
    r3r = r3;
    r3r(r3<0) = inf;
    dr3 = min(r3r,[],2);
    dr3(linear_inds) = r3lin;
    
    % If there are any inf values left, just use the closest
    % point to 0, monotonicity may have been broken by the
    % regression or loss of precision somewhere
    infv = isinf(dr3);
    if any(infv)                    
        minv = min(abs(r3(infv,:)),[],2);
        dr3(infv) = minv;
    end
    
    dr3 = reshape(dr3,size(indices));

    % Get inverse values by adding back the breakpoint offset
    def_interp = dr3 + pp.breaks(indices);
    
    def_interp = permute(def_interp,[2 3 1]);
end

boundary_points_flat = def_interp;

flat_boundary_inds = N_csum(2:end);


function [x,linear_inds] = cs_inv2(p)
% Find the roots of either linear or cubic polynomials to invert the
% interpolation
% 
% Polynomials take the form y = ax^3 + bx^2 + cx + d 

% We only care about positive real roots, so use -1 to flag either that
% there is no root (linear), or that it is complex (cubic)
x = -ones(size(p,1),3);

% Coefficients
a = p(:,1);
b = p(:,2);
c = p(:,3);
d = p(:,4);

% First get linear polynomial data
linear_inds = (a == 0) & (b == 0);

% To avoid numerical instability, remove small leading coefficients
quadratic_inds = abs(a) < 1e-9 & (b ~= 0);
cubic_inds = ~(linear_inds | quadratic_inds);

if any(linear_inds)
    x(linear_inds,1) = -d(linear_inds)./c(linear_inds);
end

if any(quadratic_inds)
    q = -1/2*(c + sign(c).*sqrt(c.^2 - 4*b.*d));
    x(quadratic_inds,1) = q(quadratic_inds)./b(quadratic_inds);
    x(quadratic_inds,2) = d(quadratic_inds)./q(quadratic_inds);
end

if any(cubic_inds)
    % Make leading coefficient equal to 1 (same solution, simpler formula) to
    % make the new formula y = x^3 + bx^2 + cx + d
    a = a(cubic_inds);
    b = b(cubic_inds)./a;
    c = c(cubic_inds)./a;
    d = d(cubic_inds)./a;

    x1 = x(cubic_inds,1); 
    x2 = x1; % empty so copy
    x3 = x1;

    % The following algorithm is from Chapter 5.6 of Numerical Recipes
    % The book notes that the equations are set up to mimimize round off error
    Q = (b.^2 - 3*c)/9;
    R = (2*b.*b.*b - 9*b.*c + 27*d)/54;

    Q3 = Q.*Q.*Q; % a lot faster than Q.^3...
    R2 = R.^2;

    % Three real roots case
    r3 = R2 < Q3;

    sQ = -2*sqrt(Q(r3));
    theta = acos(R(r3)./sqrt(Q3(r3)));
    b3 = b(r3)/3;

    x1(r3) = sQ.*cos(theta/3) - b3;
    x2(r3) = sQ.*cos((theta+2*pi)/3) - b3;
    x3(r3) = sQ.*cos((theta-2*pi)/3) - b3;

    % One real root and two complex conjugates case
    if any(~r3)
        R = R(~r3);
        b = b(~r3);

        A = -sign(R).*(abs(R) + sqrt(R2(~r3) - Q3(~r3))).^(1/3);

        B = Q(~r3)./A;
        B(A==0) = 0;

        x1(~r3) = (A + B) - b/3;

        % Imaginary part of roots is zero giving one more real root
        % (Not sure that this is possible/necessary?)
        r2 = A == B;
        if any(r2)
            x2r = x2(~r3);
            x2r(r2) = -(A(r2) + B(r2))/2 - b(r2)/3;
            x2(~r3) = x2r;
        end
    end

    x(cubic_inds,1) = x1;
    x(cubic_inds,2) = x2;
    x(cubic_inds,3) = x3;
end