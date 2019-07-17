function [outN,imageNoise,out] = early_vision_model_NoFovea(imageBitmap,degSize,timecourse,pars,sizePx,V1Mode)
%function [imageBitmap,imageNoise] = early_vision_model(imageBitmap,degSize,timecourse,fixation,pars,sizePx,V1Mode)



if ~exist('degSize','var') || isempty(degSize)
    degSize = getDegSize;
end

if ~exist('sizePx','var') || isempty(sizePx)
    sizePx = size(imageBitmap);
end
if ~exist('timecourse','var') || isempty(timecourse)
    timecourse = 250;
end

if isscalar(timecourse)
    timecourse = ones(ceil(timecourse),1);  %ms
end

if ~exist('V1Mode','var') || isempty(V1Mode)
    V1Mode = 1;
end
if ~exist('csfSelector','var') || isempty(csfSelector)
    if V1Mode == 1
        csfSelector = 'mean';
    elseif V1Mode == 2
        csfSelector = 'standard';
    else
        csfSelector = [];
    end
end

if ~exist('pars','var') || isempty(pars)
    pars = getPars(length(timecourse),V1Mode);
end

noiseConst = pars(1);
noiseFactor= pars(2);
pars(8:9) = round(pars(8:9));


if any(pars<0) 
    error('You passed negative parameters to the model. This is not allowed!');
end

imageBitmap = watson_filter_eye(imageBitmap,degSize,4);

[imageBitmap,degSize] = foveate_image0(imageBitmap,degSize,timecourse,csfSelector);

imageBitmap = imresize(imageBitmap,sizePx);

[outN,~,out] = V1(imageBitmap,degSize,4,pars(3:14),[]);

imageNoise = noiseConst+noiseFactor.*abs(out);
imageNoise = imageNoise./prod(degSize)*numel(imageBitmap).*pars(8)./12.*pars(9)./8;

