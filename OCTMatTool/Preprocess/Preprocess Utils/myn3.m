function [biasField,v_corrected,spMat] = myn3(v,u,mask,header,params,spMat)
% Inhomogeneity correction of 2D OCT data
%
% Algorithm uses a multiplicative model of inhomogeneity v = u*f where v is
% the input data, u is the corrected image, and f is the bias field.
% Algorithm operates on log transformed data so we have v~ = u~ + f~ where
% v~ = log(v)
%
% Note: currently works in 2D only
%
% Inputs:
%   v - volume to apply inhomogeneity correction to
%   u - initial estimate of corrected data
%   mask - mask for the data to run the correction on (data outside the
%          masked region is ignored)
%   header - header of the OCT data which must include pixel sizes
%   params - algorithm parameters - a struct with the fields
%               numIterations - number of iterations to run the algorithm
%                   (default: 50)
%               controlPointDistance - control point spacing in units of
%                   micron for the B-spline smoother; can be a vector
%                   [sy,sx] providing the spacing in the y and x
%                   directions, respectively (default: 200)
%               lambda - regularization parameter for B-spline fit; note
%                   that regularization can be enforced by both increasing
%                   controlPointDistance and by increasing lambda. The y
%                   and x directions can be regularized separately by input
%                   of a 2x1 vector of lambda values for the y and x
%                   directions, respectively. (default: 1e-7)
%               numHistogramBins - number of histogram bins for the
%                   histogram (default: 200)
%               convergenceThresh - threshold for convergence using the
%                   standard deviation of the difference in field estimates
%                   (default: 0.001)
%               fwhm - full width half maximum of the incremental Gaussian
%                   distribution of the bias field (default: 0.15)
%               filterNoise - noise value for the Wiener filter (default:
%                   0.015)
%
% Code is based on the N3 implementation from the authors of the N4ITK
% method - code available at:
% http://www.insight-journal.org/browse/publication/640
%
% Code for the B-spline smoothing is based on the work:
%   P.H.C. Eilers, I.D. Currie, M. Durban. "Fast and compact smoothing on
%   large multidimensional grids." Computational Statistics & Data
%   Analysis, vol. 50, pp. 61-76, 2006.

% TODO:
%   Fix input handling of parameters
%   Extend to 3D

if nargin < 6
    spMat = [];
end
if nargin < 5
    params.numIterations = 25;
    params.numHistogramBins = 200;
    params.convergenceThresh = 0.00001;
    params.fwhm = 0.15;
    params.filterNoise = 0.01;
    params.controlPointDistance = 155;
    params.lambda = 1e-7;
    params.normalize_field = false;
    params.resize_width = 256;
end

m_MaximumNumberOfIterations = params.numIterations;
m_NumberOfHistogramBins = params.numHistogramBins;
m_ConvergenceThreshold = params.convergenceThresh;
m_BiasFieldFullWidthAtHalfMaximum = params.fwhm;
m_WeinerFilterNoise = params.filterNoise;
m_controlPointDistance = params.controlPointDistance;
lambda = params.lambda;
normalize_field = params.normalize_field;
resize_width = params.resize_width;

if length(lambda) == 1
    lambda = [lambda lambda]; % [lambda_y lambda_x]
end

v = double(v);
u = double(u);

% Initial downsample
full_size = size(v);
new_size = [full_size(1),resize_width];
v_rs = imresize(v,new_size);
u_rs = imresize(u,new_size);
header.ScaleX = header.ScaleX*full_size(2)/new_size(2);

if isempty(mask)
    mask = true(size(v_rs));
else
    mask = imresize(mask,new_size);
end

% Control point distance in pixels
res = [(header.ScaleZ*1000) (header.ScaleX*1000)];
m_controlPointDistance = m_controlPointDistance./res;

% Calculate the log of the input image.
v_rs(v_rs<0) = 0;
u_rs(u_rs<0) = 0;
logInputImage = log(v_rs);
logInputU = log(u_rs);

% Remove possible nans/infs from the log input image.
logInputImage(isnan(logInputImage) | isinf(logInputImage)) = 0;
logInputU(isnan(logInputU) | isinf(logInputU)) = 0;

% Provide an initial log bias field of zeros
logBiasField = logInputImage-logInputU;
% logBiasField = zeros(size(logInputImage));

if isempty(spMat)
    % Initialize smoothing matrices
    spMat = eilersSplineBuildMatrices(m_controlPointDistance.*res,size(logInputImage),lambda,res,mask);
end

% Update spline matrices with mask information
spMat = eilersSplineUpdateMatricesNewMask(spMat,mask);

% Iterate until convergence or iterative exhaustion.
isConverged = false;
m_ElapsedIterations = 0;
while(~isConverged && m_ElapsedIterations < m_MaximumNumberOfIterations)    

    % Estimate u as v - f    
    subt = logInputImage - logBiasField;
    
    % Estimate E[u|v] after sharpening
    logSharpenedImage = sharpenImage(subt,mask,m_NumberOfHistogramBins,...
        m_BiasFieldFullWidthAtHalfMaximum,m_WeinerFilterNoise);
    
    % Estimate f
    subt2 = logInputImage - logSharpenedImage;
    
    % Smooth estimate of f  
    newLogBiasField = eilersSplineSolve(spMat,subt2);
    
    % Convergence criteria
    m_CurrentConvergenceMeasurement =...
        calculateConvergenceMeasurement(logBiasField,newLogBiasField,mask);    
    isConverged = m_CurrentConvergenceMeasurement < m_ConvergenceThreshold;

    logBiasField = newLogBiasField;
    
    m_ElapsedIterations = m_ElapsedIterations + 1;
end
if params.verbose
    fprintf('Converged in %d iterations with convergence value %2.5f\n',...
        m_ElapsedIterations,m_CurrentConvergenceMeasurement);
end

% Undo log transform
biasField = exp(logBiasField);

if normalize_field
    % Normalize field to have mean of one inside mask
    mv = mean(biasField(mask));
    biasField = biasField/mv;
end

% if test
%     % Rescale so the average intensity of the image doesn't change
%     mean_in = mean(u_rs(:));
%     mean_out = mean(v_rs(:)./biasField(:));
%     biasField = biasField*mean_out/mean_in;
% end

biasField = imresize(biasField,full_size);

v_corrected = v./biasField;


function sharpenedImage = sharpenImage(unsharpenedImage,mask,...
    m_NumberOfHistogramBins,m_BiasFieldFullWidthAtHalfMaximum,...
    m_WeinerFilterNoise)

% Build the histogram for the uncorrected image. Note that variables
% in real space are denoted by a single uppercase letter whereas their
% frequency counterparts are indicated by a trailing lowercase 'f'.
binMaximum = max(unsharpenedImage(mask));
binMinimum = min(unsharpenedImage(mask));

histogramSlope = (binMaximum - binMinimum)/(m_NumberOfHistogramBins-1);

% Create the intensity profile (within the masked region, if applicable)
% using a triangular parzen windowing scheme.
cidx = (unsharpenedImage-binMinimum)/histogramSlope + 1;
% cidx = (unsharpenedImage(mask)-binMinimum)/histogramSlope + 1;
cidxm = cidx(mask);
idx = floor(cidxm);
offset = cidxm - idx;

H1 = accumarray(idx+1,offset);
H2 = accumarray(idx,1-offset);
if max(idx) < m_NumberOfHistogramBins
    d = m_NumberOfHistogramBins-max(idx);
    H1 = cat(1,H1,zeros(d,1));
    H2 = cat(1,H2,zeros(d,1));
end
H = H1(1:end-1)+H2;
if length(H) ~= m_NumberOfHistogramBins
    error('Histogram creation failed!')
end

% % Test loop to make sure H is correct
% H = zeros(m_NumberOfHistogramBins,1);
% for i = 1:(m_NumberOfHistogramBins-1)
%     idx_i = idx == i;
%     H(i) = H(i) + sum(1-offset(idx_i));
%     H(i+1) = sum(offset(idx_i));
% end
% H(end) = H(end) + sum(idx==m_NumberOfHistogramBins);

% Determine information about the intensity histogram and zero-pad
% histogram to a power of 2.
exponent = ceil(log(m_NumberOfHistogramBins)/log(2)) + 1;
paddedHistogramSize = 2^exponent;
histogramOffset = floor((paddedHistogramSize-m_NumberOfHistogramBins)/2);
V = zeros(paddedHistogramSize,1);
V((1:m_NumberOfHistogramBins)+histogramOffset) = H;
Vf = fft(V);

% Create the Gaussian filter.
scaledFWHM = m_BiasFieldFullWidthAtHalfMaximum/histogramSlope;
expFactor = 4*log(2)/scaledFWHM^2;
scaleFactor = 2*sqrt(log(2)/pi)/scaledFWHM;
F = zeros(paddedHistogramSize,1);
F(1) = scaleFactor;
halfSize = floor(paddedHistogramSize/2);
n = 2:(halfSize+1);
F(n) = scaleFactor*exp(-(n-1).^2*expFactor);
F(paddedHistogramSize-n+2) = F(n);
if(mod(paddedHistogramSize,2) == 0)
    F(halfSize+1) = scaleFactor*exp(-0.25*paddedHistogramSize^2*expFactor);
end
Ff = fft(F,paddedHistogramSize);

% Create the Weiner deconvolution filter.
Gf = conj(Ff)./(conj(Ff).*Ff + m_WeinerFilterNoise);

% I think -> Gf should be real since F is symmetric, so why not force it to
% be real when multiplying? 
Gf = real(Gf);

Uf = Vf.*Gf;
U = ifft(Uf);
U = max(U,0);

% Compute mapping E(u|v)
n = (0:(paddedHistogramSize-1))';
numerator = (binMinimum + (n-histogramOffset)*histogramSlope).*U;
numerator = fft(numerator);
numerator = numerator.*Ff;
numerator = ifft(numerator);

denominator = U;
denominator = fft(denominator);
denominator = denominator.*Ff;
denominator = ifft(denominator);

E = real(numerator)./real(denominator);
E(isinf(E) | isnan(E)) = 0;

% Remove the zero-padding from the mapping
E = E((1:m_NumberOfHistogramBins)+histogramOffset);

% Sharpen the image with the new mapping, E(u|v)
% cidx = (unsharpenedImage-binMinimum)/histogramSlope + 1;

% sharpenedImage = interp1(1:m_NumberOfHistogramBins,E,cidx,'linear','extrap');

F = griddedInterpolant(1:m_NumberOfHistogramBins,E,'linear');
F.ExtrapolationMethod = 'linear';
sharpenedImage = reshape(F(cidx(:)),size(cidx));


% % function spMat = eilersSplineBuildMatrices2(m_controlPointDistance,fieldSize,lambda,res,mask)
% % % Pre-allocate matrices that can be reused at each iteration of the
% % % algorithm - see Eilers et al. (2006) for details
% % 
% % d = round(m_controlPointDistance);
% % 
% % % Data points
% % Y = 1:fieldSize(1);
% % X = 1:fieldSize(2);
% % 
% % % Place spline control points over the data
% % 
% % % Spline control points centered on the data
% % sp1 = round(mod(fieldSize(1),d(1)*2)/2);
% % Yc = sp1:d(1):fieldSize(1);
% % sp2 = round(mod(fieldSize(2),d(2)*2)/2);
% % Xc = sp2:d(2):fieldSize(2);
% % 
% % % Cubic B-splines are actually centered on the 3rd control point, so add
% % % 2 control points on each end to have a B-spline basis function centered
% % % on each of the above control points. Alternatively you can use repeated
% % % control points, but I think that makes things less smooth.
% % % Yc = cat(2,Yc(1)-[2*d(1) d(1)],Yc,Yc(end)+[d(1) 2*d(1)]);
% % % Xc = cat(2,Xc(1)-[2*d(2) d(2)],Xc,Xc(end)+[d(2) 2*d(2)]);
% % Yc = cat(2,Yc(1)-[3*d(1) 2*d(1) d(1)],Yc,Yc(end)+[d(1) 2*d(1) 3*d(1)]);
% % Xc = cat(2,Xc(1)-[3*d(2) 2*d(2) d(2)],Xc,Xc(end)+[d(2) 2*d(2) 3*d(2)]);
% % 
% % % N3 knot placement (see TBSpline.cc)
% % n = ceil((fieldSize-1)./(d*(1-eps)))+3;
% % start = 0.5*((fieldSize-1)-d.*(n+3));
% % Yc = (round(start(1)) + d(1)*(0:(n(1)+3)))+1;
% % Xc = (round(start(2)) + d(2)*(0:(n(2)+3)))+1;
% % 
% % % B-spline matrix in X-direction
% % order = 3;
% % knots = Xc;
% % sp_x = fastBSpline(knots,ones(length(knots)-order-1,1));
% % B2 = sp_x.getBasis(X);
% % 
% % % Keep these in temporarily for checks before removing fastBSpline
% % % dependency
% % B2b = bspl(X',Xc(3),Xc(end-2),length(Xc)-5,3);
% % B2b(:,1) = []; B2b(:,end) = [];
% % 
% % % B-spline matrix in Y-direction
% % knots = Yc;
% % sp_y = fastBSpline(knots,ones(length(knots)-order-1,1));
% % B1 = sp_y.getBasis(Y);
% % 
% % B1b = bspl(Y',Yc(3),Yc(end-2),length(Yc)-5,3);
% % B1b(:,1) = []; B1b(:,end) = [];
% % 
% % if (max(abs(B1(:)-B1b(:))) > 1e-7) || (max(abs(B2(:)-B2b(:))) > 1e-7)
% %     error('bspl inconsistent with fastBSpline')
% % end
% % 
% % [m1,n1] = size(B1);
% % [m2,n2] = size(B2);
% % 
% % % Penalty matrix
% % E1 = eye(n2); E2 = eye(n1);
% % D1 = diff(E1,2); D2 = diff(E2,2);
% % D11 = diff(E1,1); D22 = diff(E2,1);
% % 
% % % Replicate boundary conditions (does this matter? First and last B-splines
% % % are outside of the image)
% % D1 = cat(1,[-1 1 zeros(1,n2-2)], D1, [zeros(1,n2-2) 1 -1]);
% % D2 = cat(1,[-1 1 zeros(1,n1-2)], D2, [zeros(1,n1-2) 1 -1]);
% % D11 = cat(1,zeros(1,n2), D11);
% % D22 = cat(1,zeros(1,n1), D22);
% % 
% % P1 = kron(D1'*D1,E2);
% % P2 = kron(E1,D2'*D2);
% % P3 = kron(D11'*D11,D22'*D22);
% % 
% % %%%% NOTE I DONT KNOW WHAT IS RIGHT, THE CONTROL POINTS ARE FARTHER APART
% % %%%% IN THE Y DIRECTION BECAUSE OF THE PIXEL SIZE, SO SHOULDNT THE WEIGHT
% % %%%% BE HIGHER IN THE X DIRECTION? 
% % % % Weight penalties in each direction based on anisotropy
% % % % (closer in pixels means more smoothing since the physical distance is the
% % % % same)
% % if res(1) < res(2)  % y spacing smaller than x
% %     P2 = res(2)/res(1)*P2; % weight y direction higher (more smoothing)
% % else
% %     P1 = res(1)/res(2)*P1;
% % end
% % 
% % % P1 = res(2)*P1;
% % % P2 = res(1)*P2;
% % 
% % P = lambda*(P1 + P2);
% % P = P + lambda*2*P3;
% % 
% % % Data weight matrix
% % W = ones(m1,m2);
% % W(~mask) = 0; % W = mask?
% % 
% % % Compute F matrix (described in paper)
% % F = box(B1)'*W*box(B2); % F0 in paper
% % F = reshape(F,[n1,n1,n2,n2]); % F1 - make 4-D array
% % F = permute(F,[1,3,2,4]); % F2 - permute dimensions
% % F = reshape(F,[n1*n2,n1*n2]); % F3 - back to 2-D
% % 
% % spMat.F = F;
% % spMat.P = P;
% % spMat.B1 = B1;
% % spMat.B2 = B2;
% % spMat.W = W;
% % spMat.lambda = lambda;


function spMat = eilersSplineBuildMatrices(m_controlPointDistance,fieldSize,lambda,res,mask)
% Pre-allocate matrices that can be reused at each iteration of the
% algorithm - see Eilers et al. (2006) for details

crossPen = false; % include dxdy term in addition to dx^2 dy^2 terms

% d = round(m_controlPointDistance);
d = m_controlPointDistance;

fieldSize_um = (fieldSize-1).*res;

% Data points
Y = ((1:fieldSize(1))-1)*res(1);
X = ((1:fieldSize(2))-1)*res(2);

% Place spline control points over the data

% Spline control points centered on the data
sp1 = mod(fieldSize_um(1),d(1)*2)/2;
Yc = sp1:d(1):fieldSize_um(1);
sp2 = mod(fieldSize_um(2),d(2)*2)/2;
Xc = sp2:d(2):fieldSize_um(2);

% Cubic B-splines are actually centered on the 3rd control point, so add
% 2 control points on each end to have a B-spline basis function centered
% on each of the above control points. Alternatively you can use repeated
% control points, but I think that makes things less smooth.
% Yc = cat(2,Yc(1)-[2*d(1) d(1)],Yc,Yc(end)+[d(1) 2*d(1)]);
% Xc = cat(2,Xc(1)-[2*d(2) d(2)],Xc,Xc(end)+[d(2) 2*d(2)]);
Yc = cat(2,Yc(1)-[3*d(1) 2*d(1) d(1)],Yc,Yc(end)+[d(1) 2*d(1) 3*d(1)]);
Xc = cat(2,Xc(1)-[3*d(2) 2*d(2) d(2)],Xc,Xc(end)+[d(2) 2*d(2) 3*d(2)]);

% N3 knot placement (see TBSpline.cc)
n = ceil(fieldSize_um./(d*(1-eps)))+3;
start = 0.5*(fieldSize_um-d.*(n+3));
Yc = (start(1) + d(1)*(0:(n(1)+3)));
Xc = (start(2) + d(2)*(0:(n(2)+3)));

% B-spline matrix in X-direction
B1 = bspl(X',Xc(3),Xc(end-2),length(Xc)-5,3);
B1(:,1) = []; B1(:,end) = [];

% B-spline matrix in Y-direction
B2 = bspl(Y',Yc(3),Yc(end-2),length(Yc)-5,3);
B2(:,1) = []; B2(:,end) = [];

[m1,n1] = size(B1);
[m2,n2] = size(B2);

% Second order penalty matrix in each direction
E1 = eye(n1);
D1 = diff(E1,2); % x-dir
E2 = eye(n2);
D2 = diff(E2,2); % y-dir

% Replicate boundary conditions
D1 = cat(1,[-1 1 zeros(1,n1-2)], D1, [zeros(1,n1-2) 1 -1]);
D2 = cat(1,[-1 1 zeros(1,n2-2)], D2, [zeros(1,n2-2) 1 -1]);

% % Second order derivative -> divide by d^2
% % We are taking the derivative of coefficients as opposed to the spline
% % itself, so maybe this doesn't make sense
% D1 = D1/d(1)^2;
% D2 = D2/d(2)^2;

% Final penalty term
P1 = kron(D1'*D1,E2);
P2 = kron(E1,D2'*D2);
P = lambda(2)*P1 + lambda(1)*P2; % lambda = [lambda_y lambda_x]

if crossPen
    % First order penalty matrix
    D11 = diff(E1,1); 
    D22 = diff(E2,1);
    % Replicate boundary conditions
    D11 = cat(1,zeros(1,n1), D11);
    D22 = cat(1,zeros(1,n2), D22);
%     % First order derivative -> divide by d
%     D11 = D11/d(1);
%     D22 = D22/d(2);
    
    P3 = kron(D11'*D11,D22'*D22);
    
    P = P + lambda(3)*2*P3;
end

spMat.boxB1 = box(B1);
spMat.boxB2 = box(B2);
spMat.P = P;
spMat.B2 = B1;
spMat.B1 = B2;
spMat.lambda = lambda;
spMat.res = res;


function spMat = eilersSplineUpdateMatricesNewMask(spMat,mask)

[m1,n1] = size(spMat.B2);
[m2,n2] = size(spMat.B1);

% Data weight matrix
W = ones(m2,m1);
W(~mask) = 0;

% Compute F matrix (described in paper)
F = spMat.boxB2'*W*spMat.boxB1; % F0 in paper
F = reshape(F,[n2,n2,n1,n1]); % F1 - make 4-D array
F = permute(F,[1,3,2,4]); % F2 - permute dimensions
F = reshape(F,[n1*n2,n1*n2]); % F3 - back to 2-D

% Build final design matrix (multiplying P by np gives an average ssd term
% and dividing by np*res_x*res_y gives the area of the smoothed surface,
% thus normalizing each term)
np = nnz(W); % number of pixels
M = sparse(F + np/(np*spMat.res(1)*spMat.res(2))*spMat.P); % sparse solver is faster

spMat.M = M;
spMat.W = W;

% spMat = rmfield(spMat,{'boxB1','boxB2','P'});

function smoothedField = eilersSplineSolve(spMat,field)
% Solve for the smoothed B-spline field
%   See Eilers et al. (2006) for details
%
% Note that we make a small change compared to the paper by dividing the
% mean squared error part of the cost function by N, the number of pixels.
% Dividing the MSE by N is equivalent to multiplying lambda by N, so we do
% this when solving the system of equations: alpha = (F+lambda*N*P)^-1 * R

% F = spMat.F;
% P = spMat.P;
M = spMat.M;
B1 = spMat.B1;
B2 = spMat.B2;
W = spMat.W;

YY = field;

% n = nnz(W); % number of pixels

R = B1'*(W.*YY)*B2;
R = R(:); % Make it a vector

n1 = size(B1,2);
n2 = size(B2,2);

% Solve penalized system
% M = sparse(F + n/(n*spMat.res(1)*spMat.res(2))*P); % sparse solver is faster
% A = M\R;
A = M\R;
A = reshape(A,[n1,n2]); % Make it a matrix

smoothedField = B1*A*B2';

% % Check on penalty term (should be equal)
% pen = A(:)'*P*A(:);
% penf = spMat.lambda*sum(sum(imfilter(A,[-1 2 -1]/200^2,'replicate').^2 +...
%                             imfilter(A,[-1 2 -1]'/200^2,'replicate').^2 +...
%                             2*imfilter(A,[-1 1; 1 -1]/200^2,'replicate').^2));
%                         
% % Compare to penalty on smoothed field
% f1 = imfilter(smoothedField,[-1 2 -1]/spMat.res(2)^2,'replicate').^2;
% f2 = imfilter(smoothedField,[-1 2 -1]'/spMat.res(1)^2,'replicate').^2;
% f3 = 2*imfilter(smoothedField,[-1 1; 1 -1]/prod(spMat.res),'replicate').^2;
% penf2 = sum(f1(logical(W))+f2(logical(W))+f3(logical(W)));


function C = box(B)
% Box function described in Eilers paper
eL = ones(size(B,2),1)';
C = kron(B,eL).*kron(eL,B);


function cv = calculateConvergenceMeasurement(field1,field2,mask)
% Coefficient of variation of the field difference

field_diff = field1-field2;
field_diff = field_diff(mask);

% Undo log
field_diff = exp(field_diff);

% % Compute CV
% cv = std(field_diff)/mean(field_diff);

% Actual N3 code uses only standard deviation (mean(field_diff) converges
% to 1 since exp(0) = 1, so these are very similar)
cv = std(field_diff);


function B = bspl(x,xl,xr,ndx,deg)
% Set up the B-spline matrix for B-spline fitting. Assumes the spline
% functions are evenly spaced over the domain.
%
% Inputs:
%   x   -   points to evalute each B-spline at (column vector)
%   xl  -   center point of lower B-spline
%   xr  -   center point of upper B-spline
%   ndx -   number of intervals (number of spline points minus 1)
%   deg -   degree of spline
%
% Code extracted from the paper:
%   P.H.C Eilers and B.D. Marx. "Flexible Smoothing with B-splines and
%   Penalties." Statistical Science, vol. 11 (2), pp. 89-121, 1996.

dx = (xr-xl)/ndx;
t = xl + dx*(-deg:ndx-1); % knot points
T = (0*x+1)*t;
X = x*(0*t+1);
P = (X-T)/dx;
B = (T<=X)&(X<(T+dx));
r = [2:length(t) 1];
for k = 1:deg
    B = (P.*B+(k+1-P).*B(:,r))/k;
end

