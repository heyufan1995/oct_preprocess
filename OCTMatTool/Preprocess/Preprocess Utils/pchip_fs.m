function pp = pchip_fs(x,y,xx)
% Vectorized efficient pchip interpolation which assumes monotonic data,
% specifically implemented for flat space transformation of OCT data. Has
% about a 100x speedup over pchip for OCT volumes.
%
% See pchip for help and details
%
% Major differences other than vectorization between pchip_fs and pchip:
%   - no input checking
%   - no support for complex data
%   - input y must be monotonic
%
% Without enforcing monotonicity, the pchip code is not easily vectorizable
% due to finding points where the slope changes and enforcing a zero slope
% at those points. These points may be different for each value of y(i,:).

h = diff(x);
m = size(y,1);

% Compute slopes
del = diff(y,1,2)./repmat(h,m,1);
slopes = zeros(size(y));

% Check monotonicity
if ~(all(del(:) < 0) || all(del(:) > 0))
    error('Input to ''mypchip'' function must be monotonic!')
end

% Derivative values for shape-preserving Piecewise Cubic Hermite
% Interpolation. Inline and simplified version of d = pchipslopes(x,y,del)
% from pchip.m.
n = length(x);

%  Special case n=2, use linear interpolation.
if n==2  
    slopes = repmat(del(:,1),size(y,2));
end

k = 1:(n-2);
hs = h(k)+h(k+1);
w1 = (h(k)+hs)./(3*hs);
w2 = (hs+h(k+1))./(3*hs);

%  Slopes at interior points.
%  d(k) = weighted average of del(k-1) and del(k)
dmax = max(abs(del(:,k)),del(:,k+1));
dmin = min(abs(del(:,k)),del(:,k+1));
slopes(:,k+1) = dmin./(bsxfun(@times,w1,del(:,k)./dmax) + bsxfun(@times,w2,del(:,k+1)./dmax));

%  Slopes at end points. Set d(1) and d(n) via non-centered,
%  shape-preserving three-point formulae.
slopes(:,1) = ((2*h(1)+h(2))*del(:,1) - h(1)*del(:,2))/(h(1)+h(2));

% Check monotonicity (probably not necessary anymore but trivial
% computation)
indsc = sign(slopes(:,1)) ~= sign(del(:,1));
slopes(indsc,1) = 0;
indsc = (sign(del(:,1)) ~= sign(del(:,2))) & (abs(slopes(:,1)) > abs(3*del(:,1)));
slopes(indsc,1) = 3*del(indsc,1);

slopes(:,n) = ((2*h(n-1)+h(n-2))*del(:,n-1) - h(n-1)*del(:,n-2))/(h(n-1)+h(n-2));

% Check monotonicity (probably not necessary anymore but trivial
% computation)
indsc = sign(slopes(:,n)) ~= sign(del(:,n-1));
slopes(indsc,n) = 0;
indsc = (sign(del(:,n-1)) ~= sign(del(:,n-2))) & (abs(slopes(:,n)) > abs(3*del(:,n-1)));
slopes(indsc,n) = 3*del(indsc,n-1);

% Compute piecewise cubic Hermite interpolant to those values and slopes
pp = pwch(x,y,slopes,h,del); 
pp.dim = m;

if nargin == 3   % if values are wanted instead, provide them
   pp = ppval(pp,xx);
end
