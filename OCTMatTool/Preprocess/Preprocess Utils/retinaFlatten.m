function flatImage = retinaFlatten(img_data,shifts,interpMethod)
% Flatten by shifting A-scans up or down depending on the value in 'shifts'
% Extrapolation just adds zeros to the top or bottom of the A-scan

% Andrew Lang
% Updated to run slice by slice, which is about twice as fast for 49
% B-scans
% $Id: retinaFlatten.m,v 1.4 2015/01/12 18:52:29 andrew Exp $

[Y, X] = ndgrid(1:size(img_data,1),1:size(img_data,2));
Ym = bsxfun(@plus,Y,shiftdim(shifts,-1));

img_type = class(img_data);
img_data = double(img_data);
flatImage = zeros(size(img_data));
if strcmp(version,'7.13.0.564 (R2011b)')
    % Bug in this version for griddedInterpolant
    for i = 1:size(img_data,3)
        flatImage(:,:,i) = interp2(X,Y,double(img_data(:,:,i)),X,Ym(:,:,i),interpMethod,0);
    end
else
    for i = 1:size(img_data,3)
        F = griddedInterpolant(Y,X,img_data(:,:,i),interpMethod);
        if isprop(F,'ExtrapolationMethod')
            % For consistent behavior with the old function, set
            % extrapolation to output nans
            F.ExtrapolationMethod = 'none';
        end
        flatImage(:,:,i) = F(Ym(:,:,i),X);
    end
end
% Set extrapolatd values to 0
flatImage(isnan(flatImage)) = 0;
flatImage = cast(flatImage,img_type);

