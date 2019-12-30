function [Xr, Yr] = rot2d_center(X,Y,angle)

% Rotation matrix 
R = [cosd(angle) -sind(angle); sind(angle) cosd(angle)];

Pr = R*[X(:)'; Y(:)'];
Xr = reshape(Pr(1,:),size(X));
Yr = reshape(Pr(2,:),size(Y));