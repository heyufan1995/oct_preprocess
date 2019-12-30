function [control_pts, bd_pts] = readControlPoints(seg_file,n_ascan)
% Read control points from a mat file or OCTToolbox file. Output the
% interpolated boundaries 

% Andrew Lang
% $Id: readControlPoints.m,v 1.2 2015/01/13 14:37:31 andrew Exp $

if nargin < 2
    n_ascan = [];
end

[pathstr,name,ext] = fileparts(seg_file);

if strcmp(ext,'.txt')
    % OCT Toolbox file
    control_pts = readOTBControlPoints(seg_file,[],false);
    
    for i = 1:size(control_pts,1)
        for j = 1:size(control_pts,2)
            % OCT toolbox indexing starts at 0
            control_pts{i,j} = control_pts{i,j} + 1;
        end
    end    
elseif strcmp(ext,'.mat')
    % Saved from octViewer
    S = load(seg_file);
    
    if isfield(S,'bd_pts')
        bd_pts = S.bd_pts;
        control_pts = [];
        return
    else
        control_pts = S.control_pts;
    end
end

if nargout == 1
    return
end

if isempty(n_ascan)
    % Get the number of A-scans as the maximum value in the input
    ev = cellfun('isempty', control_pts);
    control_pts(ev) = {0};
    n_ascan = cellfun(@(x) max(x(:,1)),control_pts);
    n_ascan = max(n_ascan(:));
end

xpts = 1:n_ascan;
bd_pts = zeros(n_ascan,size(control_pts,1),size(control_pts,2));
for i = 1:size(control_pts,1)
    n_empty = 0; % Number of missing boundaries
    y_prev = [];
    for j = 1:size(control_pts,2)
        ctrl_pts = control_pts{i,j};

        if isequal(ctrl_pts,0) || isempty(ctrl_pts)
            if i ~= 11
                n_empty = n_empty + 1;
            end
            continue
        end
        
        % Sort for interpolate function just in case
        [~,inds] = sort(ctrl_pts(:,1));
        ctrl_pts = ctrl_pts(inds,:);

        % Interpolate control points
        pts = interpolateCtrlPts(ctrl_pts);
        x = pts(:,1);
        y = pts(:,2);
        
        % Check if we need to extrapolate the end points
        if ~isequal(xpts',x);
            y = interp1(x,y,xpts','nearest','extrap');            
        end
        
        % Make sure current boundary is equal to or below the previous
        % boundary (can happen due to the spline fit)
        if ~isempty(y_prev)
            inds = y < y_prev;
            y(inds) = y_prev(inds);
        end
        
        bd_pts(:,i,j) = y;
        y_prev = y;
    end
end

end
%%
function pts = interpolateCtrlPts(ctrl_pts)
% Interpolate a vector of control points to get values between the first and
% last point
    
    % Assume sorted
    xmin = ctrl_pts(1,1);
    xmax = ctrl_pts(end,1);
    
    xpts = (xmin:xmax)';
    
    % remove duplicate points 
    [pt1,m,n] = unique(ctrl_pts(:,1));
    pt2 = ctrl_pts(m,2);
    ctrl_pts = [pt1, pt2];
    
%     pts = interp1(ctrl_pts(:,1),ctrl_pts(:,2),xpts,'cubic',nan);
    pts = cubicSplineInterp(ctrl_pts(:,1),ctrl_pts(:,2),xpts);
    pts = cat(2,xpts,pts);
end

%% 
function splPts = cubicSplineInterp(xpts,ypts,xptsOut)

    xpts = xpts-xpts(1)+1;
    n = length(xpts);
    
    Y = calculateNaturalCubicSpline(n-1,ypts);

    splPts = zeros(size(xptsOut));
    for i = 1:n-1
        splPts(xpts(i)) = ypts(i);

        t = (((xpts(i)+1):(xpts(i+1)-1))-xpts(i))/(xpts(i+1)-xpts(i));
        splPts((xpts(i)+1):(xpts(i+1)-1)) = Y(i,1) + Y(i,2)*t + Y(i,3)*t.*t + Y(i,4)*t.*t.*t;
    end
    splPts(xpts(n)) = ypts(n);

end

function C = calculateNaturalCubicSpline(n,x)

    gamma = zeros(n+1,1);
    delta = zeros(n+1,1);
    D = zeros(n+1,1);

    gamma(1) = 0.5;
    for i = 2:n
        gamma(i) = 1/(4-gamma(i-1));
    end
    gamma(end) = 1/(2-gamma(end-1));

    delta(1) = 3*(x(2)-x(1))*gamma(1);
    for i = 2:n
        delta(i) = (3*(x(i+1)-x(i-1))-delta(i-1))*gamma(i); 
    end
    delta(end) = (3*(x(end)-x(end-1))-delta(end-1))*gamma(end);

    D(end) = delta(end);
    for i = n:-1:1
        D(i) = delta(i) - gamma(i)*D(i+1);
    end

    C = zeros(n,4);
    for i = 1:n
        C(i,:) = [x(i),D(i),3*(x(i+1)-x(i))-2*D(i)-D(i+1),2*(x(i)-x(i+1))+D(i)+D(i+1)];
    end
    
end