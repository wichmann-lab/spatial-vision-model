function imageBitmap = watson_filter_eye(imageBitmap, degSize, pupil)
% function imageBitmap = watson_filter_eye(imageBitmap, deg_size, pupil)
% this function computes an image filtered by the eyes MTF for an optimally
% corrected young observer
% Pupil is the pupil diameter in mm.

global useGPU
if isempty(useGPU)
    useGPU = false;
end

persistent MTFs
persistent degSizes
persistent imSizes
persistent pupils

if ~exist('pupil','var') || isempty(pupil)
    pupil = watson_pupilSize(mean(imageBitmap(:)).*prod(degSize),25);%mm diameter
end

%% image padding (DEACTIVATED FOR SPEED)
% imageBitmapPad = pad_image(imageBitmap,1);
% imSizePad = size(imageBitmapPad);
imSize = size(imageBitmap);
imSize = imSize(1:2);


%% search for saved values
found = false;
if ~isempty(MTFs)
    for iMTF = 1:length(MTFs)
        if all(degSize==degSizes{iMTF}) && (numel(size(imageBitmap)) == numel(imSizes{iMTF})) && all(size(imageBitmap)==imSizes{iMTF}) && (pupil == pupils{iMTF})
            found = true;
            idxMTF = iMTF;
        end
    end
end

if found
    MTF = MTFs{idxMTF};
else
    idxMTF = length(MTFs)+1;
    % calculate frequencies
    %% calculate filters
    
    x = 1:imSize(2);  % for frequency space
    y = 1:imSize(1);
    x = (x-ceil(mean(x)))./degSize(2);  % normalize to mean in middle, and freq in cyc/deg
    y = (y-ceil(mean(y)))./degSize(1);
    
    freq = sqrt(bsxfun(@plus ,x.^2,(y').^2));
    MTF  = watson_MTF(freq,pupil);
    MTF  =ifftshift(MTF);
    if useGPU 
        MTF = gpuArray(MTF);
    end
    MTFs{idxMTF} = MTF;
    degSizes{idxMTF} = degSize;
    imSizes{idxMTF} = size(imageBitmap);
    pupils{idxMTF} = pupil;
end

% compute fourier transform
imageFourier = fft2(imageBitmap);
imageFourier = imageFourier .*MTF;

imageBitmap  = ifft2(imageFourier,'symmetric'); % seems to work so far... I hope this is no error...

imageBitmap = imageBitmap.*pupil.^2.*pi./4; % conversion to Troland

