function [biasField,v_corrected] = myn3_3d(v,u,mask,header,params)
% Full 3D implementation of N3, requires a lot of memory to run

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
    lambda = [lambda lambda lambda]; % [lambda_y lambda_x lambda_z]
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
res = [(header.ScaleZ*1000) (header.ScaleX*1000) (header.Distance*1000)];
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

% Initialize smoothing matrices
spMat = eilersSplineBuildMatrices(m_controlPointDistance.*res,size(logInputImage),lambda,res,mask);

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
fprintf('Converged in %d iterations with convergence value %2.5f\n',...
    m_ElapsedIterations,m_CurrentConvergenceMeasurement);

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

biasField = imresize(biasField,full_size(1:2));

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
cidx = (unsharpenedImage(mask)-binMinimum)/histogramSlope + 1;
idx = floor(cidx);
offset = cidx - idx;

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
cidx = (unsharpenedImage-binMinimum)/histogramSlope + 1;
sharpenedImage = interp1(1:m_NumberOfHistogramBins,E,cidx,'linear','extrap');


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
Z = ((1:fieldSize(3))-1)*res(3);

% Place spline control points over the data
% N3 knot placement (see TBSpline.cc)
n = ceil(fieldSize_um./(d*(1-eps)))+3;
start = 0.5*(fieldSize_um-d.*(n+3));
Yc = (start(1) + d(1)*(0:(n(1)+3)));
Xc = (start(2) + d(2)*(0:(n(2)+3)));
Zc = (start(3) + d(3)*(0:(n(3)+3)));

% B-spline matrix in Y-direction
B1 = bspl(Y',Yc(3),Yc(end-2),length(Yc)-5,3);
B1(:,1) = []; B1(:,end) = [];

% B-spline matrix in X-direction
B2 = bspl(X',Xc(3),Xc(end-2),length(Xc)-5,3);
B2(:,1) = []; B2(:,end) = [];

% B-spline matrix in Z-direction
B3 = bspl(Z',Zc(3),Zc(end-2),length(Zc)-5,3);
B3(:,1) = []; B3(:,end) = [];

[m1,n1] = size(B1);
[m2,n2] = size(B2);
[m3,n3] = size(B3);

% Second order penalty matrix in each direction
E1 = eye(n1);
D1 = diff(E1,2); % y-dir
E2 = eye(n2);
D2 = diff(E2,2); % x-dir
E3 = eye(n3);
D3 = diff(E3,2); % z-dir

% Replicate boundary conditions
D1 = cat(1,[-1 1 zeros(1,n1-2)], D1, [zeros(1,n1-2) 1 -1]);
D2 = cat(1,[-1 1 zeros(1,n2-2)], D2, [zeros(1,n2-2) 1 -1]);
D3 = cat(1,[-1 1 zeros(1,n3-2)], D3, [zeros(1,n3-2) 1 -1]);

% Final penalty term (lambda = [lambda_y lambda_x lambda_z])
P1 = kron(kron(D1'*D1,E2),E3);
P2 = kron(kron(E1,D2'*D2),E3);
P3 = kron(kron(E1,E2),D3'*D3);
P = lambda(1)*P1 + lambda(2)*P2 + lambda(3)*P3; 

% if crossPen
%     % First order penalty matrix
%     D11 = diff(E1,1); 
%     D22 = diff(E2,1);
%     % Replicate boundary conditions
%     D11 = cat(1,zeros(1,n1), D11);
%     D22 = cat(1,zeros(1,n2), D22);
% %     % First order derivative -> divide by d
% %     D11 = D11/d(1);
% %     D22 = D22/d(2);
%     
%     P3 = kron(D11'*D11,D22'*D22);
%     
%     P = P + lambda(3)*2*P3;
% end

% Data weight matrix
W = ones(m1,m2,m3);
W(~mask) = 0;

% Compute F matrix (described in paper)
F = rho(box(B1),W,1);
F = rho(box(B2),F,2);
F = rho(box(B3),F,3);
F = reshape(F,[n1,n1,n2,n2,n3,n3]);
F = permute(F,[1 3 5 2 4 6]);
F = reshape(F,[n1*n2*n3,n1*n2*n3]);

% Build final design matrix (multiplying P by np gives an average ssd term
% and dividing by np*res_x*res_y gives the area of the smoothed surface,
% thus normalizing each term)
np = nnz(W); % number of pixels
M = sparse(F + np/(np*prod(res))*P); % sparse solver is faster

spMat.M = M;
spMat.P = P;
spMat.B1 = B1;
spMat.B2 = B2;
spMat.B3 = B3;
spMat.W = W;
spMat.lambda = lambda;
spMat.res = res;


function smoothedField = eilersSplineSolve(spMat,field)
% Solve for the smoothed B-spline field
%   See Eilers et al. (2006) for details
%
% Note that we make a small change compared to the paper by dividing the
% mean squared error part of the cost function by N, the number of pixels.
% Dividing the MSE by N is equivalent to multiplying lambda by N, so we do
% this when solving the system of equations: alpha = (F+lambda*N*P)^-1 * R

M = spMat.M;
B1 = spMat.B1;
B2 = spMat.B2;
B3 = spMat.B3;
W = spMat.W;

YY = field;

R = rho(B1,YY.*W,1);
R = rho(B2,R,2);
R = rho(B3,R,3);

n1 = size(B1,2);
n2 = size(B2,2);
n3 = size(B3,2);

% Solve penalized system
A = M\R(:);
A = reshape(A,[n1,n2,n3]); % Make it a matrix

smoothedField = rho(B3',rho(B2',rho(B1',A,1),2),3);


function C = rho(A,B,p)
% Rho function described in Eilers paper for 3d fit
sa = size(A);
sb = size(B);
n = length(sb);
ip = [(p+1):n, 1:(p-1)];
B = permute(B,[p,ip]);
B = reshape(B,[sb(p),prod(sb(ip))]);
C = A'*B;
C = reshape(C,[sa(2),sb(ip)]);
C = ipermute(C,[p,ip]);


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

