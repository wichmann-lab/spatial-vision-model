function [imageBitmap,degSize] = cutFovea(imageBitmap,degSize,fixation,sizePx,refLum)
%function imageBitmap = cutFovea(imageBitmap,degSize,fixation,sizePx,refLum)
%
% This is a functional form of the first step in models without a proper
% retina.

global useGPU
if isempty(useGPU)
    useGPU = false;
end
persistent window % The foveal window. This is given -> never changes
persistent savedSize

degSizeFinal = [2,2]; % according to Anita Hendrickson (2005) "Organization of the adult primate fovea" in "Macular degeneration"
TukeyFactor = 1; % foveola region really flat == 1 ° with cosine this should hit the 1°20' of the chapter.


%% Cut out Fovea/fill with 0s
if any(degSizeFinal> degSize)
    imageBitmap = padarray(imageBitmap,round(0.5*max(0,(degSizeFinal-degSize).*size(imageBitmap)./degSize)));
    degSize(degSize<degSizeFinal) = degSizeFinal(degSize<degSizeFinal);
end
if any(degSizeFinal< degSize)
    x = linspace(0,degSize(2),size(imageBitmap,2))-fixation(2);
    y = linspace(0,degSize(1),size(imageBitmap,1))'-fixation(1);
    imageBitmap = imageBitmap(abs(x)<degSizeFinal(1)/2,abs(y)<degSizeFinal(1)/2);
    degSize = degSizeFinal;
end

if useGPU
    if numel(sizePx)==1
        imResized = imresize(imageBitmap,sizePx./size(imageBitmap,2));
    else % will break in too old MATLAB
        imResized = imresize(imageBitmap,sizePx);
    end
else
    if numel(sizePx)==1
        imResized = imresize(imageBitmap,[sizePx,sizePx]);
    else
        imResized = imresize(imageBitmap,sizePx);
    end
end

if isempty(window) || ~all(savedSize==sizePx)
    window = image_window_tukey(size(imResized),degSizeFinal,degSizeFinal/2,degSizeFinal/2,0,TukeyFactor);
    if useGPU
        window = gpuArray(window);
    end
    savedSize = sizePx;
end

if ~exist('refLum','var') || isempty(refLum)
    imageBitmapW = imResized .* window;
    refLum = sum(imageBitmapW(:))./sum(window(:));
end

imageBitmap = ((imResized./refLum)-1).*window;

