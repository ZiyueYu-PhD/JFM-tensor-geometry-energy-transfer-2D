function [y,H]=sharpfilt2(x,Wn,sigma,lowpass,padding)
% Usage: [y,H]=sharpfilt2(x,Wn,sigma,lowpass,padding)
% Given a two-dimensional array "x", gaussfilt2 filters array "x",
% returning the results in "y". The filter is constructed with a 
% nominal cutoff frequency "Wn" (normalized by the Nyquist frequency = half 
% the sampling frequency) and conditioned using a Gaussian window with 
% width "sigma" (normalized by the array width). If lowpass=0, a high-pass
% filter is constructed; otherwise it's low-pass. Both Wn and sigma can 
% affect the effective cutoff frequency; make sure sigma is big enough. The 
% filter's frequency transfer function is returned in "H".

% Written  April 2022 by Xinyu Si 

% Input variables:
% x: two dimensional array
% Wn: nominal cutoff frequency (normalized by Nyquist frequency)
% sigma: Gaussian window width
% lowpass: default to be 1, if 0, construct high-pass filter
% padding: padding method, choose from a constant, 'symmetric',
% 'replicate', and 'circular'

% Output values
% y: filterd x
% H: filter'f frequency transfer function

lowpassdefault=1;
Htoosmall=0.95; % warn if transfer fcn peak isn't at least this high

if nargin<2
    error(['Usage: [y,H] = ' mfilename '(x,Wn,sigma,[lowpass])'])
end

if ~exist('lowpass','var') || isempty(lowpass)
    lowpass=lowpassdefault;
end

% -=- Build filter -=-
if Wn>1
    error('MATLAB:gaussfilt3:freqExceedsNyquist', ...
        'Frequency Wn exceeds Nyquist limit for given array dimensions.')
end

% build 3D freq space
[w1,w2]=freqspace(size(x,[1,2])); 
[W1,W2] = meshgrid(w1,w2);
freqamp=sqrt(W1.^2 + W2.^2);

Hd=ones(size(x));

if lowpass
    Hd(freqamp>Wn)=0; % enforce cutoff frequency
else % highpass!
    Hd(freqamp<Wn)=0; % enforce cutoff frequency
end

win=exp(-(freqamp.^2)/2/sigma^2); % Gaussian window
win=win./max(win(:));

% [organize Hd to hd so that in each dimension, it goes from n = 0~N-1]
% if dimension of Hd (N) is odd -> n = 0 initially at (N+1)/2
% if dimension of Hd is even -> n = 0 initially at N/2 + 1
% hd = flip(rot90(Hd,2)); % equivalent to reshape(flipud(Hd(:)),size(Hd))
hd = reshape(flipud(Hd(:)),size(Hd));% n = 0 now at position (N+1)/2 [odd] N/2 - 1 [even]
hd = fftshift(hd); % n = 0 now at position N for both
hd = reshape(flipud(hd(:)),size(hd)); % n = 0 now at position 1 for both

% [convert hd in to spacial domain]
h = fftshift(ifftn(hd));
h = reshape(flipud(h(:)),size(h)); % rotate for use with imfilter

% apply the filter
h = h.*win;

% abs(H) should be the filtered Hd in frequency space
H=fftshift(fftn(h));

if (max(abs(H(:)))<Htoosmall)
    warning('MATLAB:gaussfilt3:lowAmplitudeTransferFunction', ...
        ['H peaks at a small value. ' ...
        'Increase sigma for better performance.'])
end
%% 
% % -=- Apply filter -=-
% y=imfilter(x,h); % transform back to real space
% % y=imfilter(x,h,'circular'); % transform back to real space
%%
h = reshape(flipud(h(:)),size(h)); % rotate for convolution

if rem(length(h),2)==0
    h = h(2:end,2:end);
end

hsize = size(h);
A = x;
% padding = 'circular';
outSize = size(A);

A = padImage(A, hsize, padding); % pad image
fftSize = size(A);
A = ifftn( fftn(A, fftSize) .* fftn(h, fftSize), 'symmetric' ); %perform convolution

start = 1 + size(A)-outSize;
stop = start + outSize - 1;

y = A(start(1):stop(1),start(2):stop(2));

end

function [A, padSize] = padImage(A, hsize, padding)

padSize = computePadSize(size(A), hsize);

if ischar(padding)
    method = padding;
    padVal = [];
else
    method = 'constant';
    padVal = padding;
end

A = padarray_algo(A, padSize, method, padVal, 'both');

end

function padSize = computePadSize(sizeA, sizeH)

rankA = numel(sizeA);
rankH = numel(sizeH);

sizeH = [sizeH ones(1,rankA-rankH)];

padSize = floor(sizeH/2);

end

function b = padarray_algo(a, padSize, method, padVal, direction)
%PADARRAY_ALGO Pad array.
%   B = PADARRAY_AGLO(A,PADSIZE,METHOD,PADVAL,DIRECTION) internal helper
%   function for PADARRAY, which performs no input validation.  See the
%   help for PADARRAY for the description of input arguments, class
%   support, and examples.

%   Copyright 2014 The MathWorks, Inc.

if isempty(a)
    
    numDims = numel(padSize);
    sizeB = zeros(1,numDims);
    
    for k = 1: numDims
        % treat empty matrix similar for any method
        if strcmp(direction,'both')
            sizeB(k) = size(a,k) + 2*padSize(k);
        else
            sizeB(k) = size(a,k) + padSize(k);
        end
    end
    
    b = mkconstarray(class(a), padVal, sizeB);
    
elseif strcmpi(method,'constant')
    
    % constant value padding with padVal
    b = ConstantPad(a, padSize, padVal, direction);
else
    
    % compute indices then index into input image
    aSize = size(a);
    aIdx = getPaddingIndices(aSize,padSize,method,direction);
    b = a(aIdx{:});
end

if islogical(a)
    b = logical(b);
end
end

function b = ConstantPad(a, padSize, padVal, direction)

numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
sizeB = zeros(1,numDims);
for k = 1:numDims
    M = size(a,k);
    switch direction
        case 'pre'
            idx{k}   = (1:M) + padSize(k);
            sizeB(k) = M + padSize(k);
            
        case 'post'
            idx{k}   = 1:M;
            sizeB(k) = M + padSize(k);
            
        case 'both'
            idx{k}   = (1:M) + padSize(k);
            sizeB(k) = M + 2*padSize(k);
    end
end

% Initialize output array with the padding value.  Make sure the
% output array is the same type as the input.
b         = mkconstarray(class(a), padVal, sizeB);
b(idx{:}) = a;
end

function aIdx = getPaddingIndices(aSize,padSize,method,direction)
%getPaddingIndices is used by padarray and blockproc. 
%   Computes padding indices of input image.  This is function is used to
%   handle padding of in-memory images (via padarray) as well as
%   arbitrarily large images (via blockproc).
%
%   aSize : result of size(I) where I is the image to be padded
%   padSize : padding amount in each dimension.  
%             numel(padSize) can be greater than numel(aSize)
%   method : X or a 'string' padding method
%   direction : pre, post, or both.
%
%   See the help for padarray for additional information.

% Copyright 2010 The MathWorks, Inc.

% make sure we have enough image dims for the requested padding
if numel(padSize) > numel(aSize)
    singleton_dims = numel(padSize) - numel(aSize);
    aSize = [aSize ones(1,singleton_dims)];
end

switch method
    case 'circular'
        aIdx = CircularPad(aSize, padSize, direction);
    case 'symmetric'
        aIdx = SymmetricPad(aSize, padSize, direction);
    case 'replicate' 
        aIdx = ReplicatePad(aSize, padSize, direction);
end
end

%%%
%%% CircularPad
%%%
function idx = CircularPad(aSize, padSize, direction)

numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
    M = aSize(k);
    dimNums = uint32(1:M);
    p = padSize(k);
    
    switch direction
        case 'pre'
            idx{k}   = dimNums(mod(-p:M-1, M) + 1);
            
        case 'post'
            idx{k}   = dimNums(mod(0:M+p-1, M) + 1);
            
        case 'both'
            idx{k}   = dimNums(mod(-p:M+p-1, M) + 1);
            
    end
end
end

%%%
%%% SymmetricPad
%%%
function idx = SymmetricPad(aSize, padSize, direction)

numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
    M = aSize(k);
    dimNums = uint32([1:M M:-1:1]);
    p = padSize(k);
    
    switch direction
        case 'pre'
            idx{k}   = dimNums(mod(-p:M-1, 2*M) + 1);
            
        case 'post'
            idx{k}   = dimNums(mod(0:M+p-1, 2*M) + 1);
            
        case 'both'
            idx{k}   = dimNums(mod(-p:M+p-1, 2*M) + 1);
    end
end
end

%%%
%%% ReplicatePad
%%%
function idx = ReplicatePad(aSize, padSize, direction)

numDims = numel(padSize);

% Form index vectors to subsasgn input array into output array.
% Also compute the size of the output array.
idx   = cell(1,numDims);
for k = 1:numDims
    M = aSize(k);
    p = padSize(k);
    onesVector = uint32(ones(1,p));
    
    switch direction
        case 'pre'
            idx{k}   = [onesVector 1:M];
            
        case 'post'
            idx{k}   = [1:M M*onesVector];
            
        case 'both'
            idx{k}   = [onesVector 1:M M*onesVector];
    end
end
end

function out = mkconstarray(class, value, size)
%MKCONSTARRAY creates a constant array of a specified numeric class.
%   A = MKCONSTARRAY(CLASS, VALUE, SIZE) creates a constant array 
%   of value VALUE and of size SIZE.

%   Copyright 1993-2013 The MathWorks, Inc.  

out = repmat(feval(class, value), size);
end