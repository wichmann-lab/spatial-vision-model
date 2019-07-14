function [imageBitmap,degSize] = foveate_image0(imageBitmap,degSize,timecourse,csfSelector)
% function [imageBitmap,degSize] = foveate_image1(imageBitmap,degSize,timecourse,fixation,pars,sizePx,csfSelector)
% This is the foveation for the model without a proper Retina, It models
% only a normalization for the whole image and a hand tuned (spline) csf.
% This also has no foveal cut out.

global useGPU

if isempty(useGPU)
    useGPU=false;
end

if ~exist('csfSelector','var') || isempty(csfSelector)
    csfSelector = [];
end

persistent csfs
persistent timecourses
persistent savedSize
persistent savedDegSize

if ~exist('degSize','var') || isempty(degSize)
    degSize = getDegSize;
end

timecourse = timecourse(:);


%% "adaptation"

imageBitmap = imageBitmap./mean(imageBitmap(:));
imageBitmap = imageBitmap-1;

%% try to find saved csf/timecourse

found = false;
if ~isempty(csfs)
    for icsf = 1:length(csfs)
        if length(timecourses{icsf})==length(timecourse) && all(timecourse==timecourses{icsf}) && all(savedSize{icsf}==size(imageBitmap)) && all(savedDegSize{icsf}==degSize)
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
    
    
    csf = spline(log2(fgrid),sens,log2(f));
    csf(csf<0| ~isfinite(csf)) = 0;
    csf = ifftshift(csf);
    if useGPU
        csf = gpuArray(csf);
    end
    csfs{idxCsf}=csf;
    timecourses{idxCsf}=timecourse;
    savedSize{idxCsf}=size(csf);
    savedDegSize{idxCsf}=degSize;
end

%% filter with csf
imageFFT = fft2(imageBitmap);
imageFFT = imageFFT.*csf;

imageBitmap = ifft2(imageFFT,'symmetric');
