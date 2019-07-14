function [pyr,freq,pyrGradBW,orient,filterOut] = decomp_Gabor(imageBitmap,degSize,freqRange,nFreq,nOrient,bandwidth)
%function [pyr,f,orient] = decomp_Gabor(imageBitmap,degSize,freqRange,nFreq,nOrient,bandwidth)
% computes a decomposition with the given parameters with a saved
% filterbank

global useGPU

persistent filter;
persistent nOrientlocal;
persistent savedFreq;
persistent savedDegSize;
persistent savedOrient;
persistent savedbandwidth;
persistent savedImSize;
persistent filterDerivative;

if isempty(useGPU)
    useGPU = false;
end

imSize = size(imageBitmap);
freq = exp(linspace(log(freqRange(2)),log(freqRange(1)),nFreq))';
orient = linspace(0,pi,(nOrient+1));
orient = orient(1:nOrient);

%filter = filter_logGabor(size(imageBitmap),degSize,freqRange(2),orient,bandwidth);


p = exp(log(freqRange(2)/freqRange(1))/nFreq);
%[pyr,f] = pyramid(imageBitmap,filter,nFreq,p);

imageFourier = fft2(imageBitmap);
if useGPU
    pyr = nan([size(imageBitmap),nOrient,nFreq],'gpuArray');
    orient = gpuArray(orient);
    freq = gpuArray(freq);
    degSize = gpuArray(degSize);
    x = gpuArray.linspace(-imSize(2)/2/degSize(2),(imSize(2)-2)/2/degSize(2),imSize(2));
    y = gpuArray.linspace(-imSize(1)/2/degSize(1),(imSize(1)-2)/2/degSize(1),imSize(1));
    f = sqrt(bsxfun(@plus,x.^2,y'.^2));
    orient0 = bsxfun(@(x,y) atan2(y,x),x,y');
else
    pyr = nan([size(imageBitmap),nOrient,nFreq]);
    x = 1:imSize(2);
    y = 1:imSize(1);
    x = (x-ceil(mean(x)))./degSize(2);
    y = (y-ceil(mean(y)))./degSize(1);
    f = sqrt(bsxfun(@plus,x.^2,y'.^2));
    orient0 = bsxfun(@(x,y) atan2(y,x),x,y');
end

if nargout>2
    if useGPU
        pyrGradBW = nan([size(imageBitmap),nOrient,nFreq,2],'gpuArray');
    else
        pyrGradBW = nan([size(imageBitmap),nOrient,nFreq,2]);
    end
end

if ~exist('filter','var') || isempty(filter) || numel(filter)~=nFreq || isempty(nOrientlocal) || nOrientlocal ~=nOrient ...
        || isempty(savedFreq) || any(savedFreq ~=freq)|| isempty(savedDegSize) || any(savedDegSize ~=degSize)...
        || isempty(savedbandwidth) || any(savedbandwidth ~=bandwidth)|| isempty(savedOrient) || any(savedOrient ~=orient) ...
        || isempty(savedImSize) || any(savedImSize~=imSize)
    filter = cell(nFreq,1);
    nOrientlocal = nOrient;
    savedFreq = freq;
    savedDegSize = degSize;
    savedOrient = orient;
    savedbandwidth = bandwidth;
    savedImSize = imSize;
    filterDerivative = cell(nFreq,1);
end
for iFilter = 1:nFreq
    if isempty(filter{iFilter}) 
        filter{iFilter} = ifftshift(ifftshift(filter_logGabor(size(imageFourier),degSize,freq(iFilter),orient,bandwidth,f,orient0),1),2);
    end
    pyr(:,:,:,iFilter) = ifft2(bsxfun(@times,filter{iFilter},imageFourier));    
    if nargout> 2 
        if isempty(filterDerivative{iFilter}) 
            filterDerivative{iFilter} = ifftshift(ifftshift(filter_logGabor_derivative(size(imageFourier),degSize,freq(iFilter),orient,bandwidth,f,orient0),1),2);
        end
        pyrGradBW(:,:,:,iFilter,:) = ifft2(bsxfun(@times,filterDerivative{iFilter},imageFourier));
    end
    if nargout>4
        filterOut =filter;
    end
end
