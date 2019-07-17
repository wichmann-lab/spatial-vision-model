function [imageBitmap,degSize] = foveate_image1(imageBitmap,degSize,timecourse,fixation,refLum,sizePx,csfSelector)
% function [imageBitmap,degSize] = foveate_image1(imageBitmap,degSize,timecourse,fixation,refLum,sizePx,csfSelector)
% This is the foveation for the model without a proper Retina, It models
% only a normalization for the whole image and a hand tuned (spline) csf.

global useGPU

if isempty(useGPU)
    useGPU=false;
end

if ~exist('sizePx','var') || isempty(sizePx)
    sizePx = 256;
end
if ~exist('csfSelector','var') || isempty(csfSelector)
    csfSelector = [];
end
if ~exist('refLum','var') || isempty(refLum)
    refLum = [];
end
persistent csfs
persistent timecourses
persistent savedSize
persistent savedDegSize
persistent savedCsfSelector

if ~exist('degSize','var') || isempty(degSize)
    degSize = getDegSize;
end

timecourse = timecourse(:);


%% "adaptation" -> moved to cutFovea?

%imageBitmap = imageBitmap./mean(imageBitmap(:));
%imageBitmap = imageBitmap-1;


%% Cut out Fovea/fill with 0s
[imageBitmap,degSize] = cutFovea(imageBitmap,degSize,fixation,sizePx,refLum);

%% try to find saved csf/timecourse

found = false;
if ~isempty(csfs)
    for icsf = 1:length(csfs)
        if length(timecourses{icsf})==length(timecourse) && all(timecourse==timecourses{icsf}) && all(savedSize{icsf}==size(imageBitmap)) && all(savedDegSize{icsf}==degSize) && strcmp(savedCsfSelector{icsf},csfSelector)
            found = true;
            idxCsf = icsf;
        end
    end
end

if found
    csf = csfs{idxCsf};
else
    idxCsf = length(csfs)+1;
    
    
    imSize= size(imageBitmap);
    
    %% get csf
    center = ceil((imSize+1)/2);
    x = ((1:imSize(2))-center(2))./degSize(2);
    y = ((1:imSize(1))-center(1))./degSize(1);
    x = x.^2;
    y = y.^2;
    f = sqrt(bsxfun(@plus,x,y'));
    [fgrid,sens] = getCsf(csfSelector,timecourse);
    
    
    csf = interp1(log2(fgrid),sens,log2(f),'pchip');
    csf(csf<0| ~isfinite(csf)) = 0;
    csf = ifftshift(csf);
    if useGPU
        csf = gpuArray(csf);
    end
    csfs{idxCsf}=csf;
    timecourses{idxCsf}=timecourse;
    savedSize{idxCsf}=size(csf);
    savedDegSize{idxCsf}=degSize;
    savedCsfSelector{idxCsf}=csfSelector;
end

%% filter with csf
imageFFT = fft2(imageBitmap);
imageFFT = imageFFT.*csf;

imageBitmap = ifft2(imageFFT,'symmetric');
