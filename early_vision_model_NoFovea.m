function [outN,imageNoise,out] = early_vision_model3_NoFovea(imageBitmap,degSize,timecourse,pars,sizePx)
%function [imageBitmap,imageNoise] = early_vision_model3(imageBitmap,degSize,timecourse,fixation,pars,sizePx)
% This is the model in which only a few orientations enter the
% normalization pool. This should fit oblique masking additionally. 
%
% According to Itti, Koch & Brown (2000) this should also fit orientation
% and frequency discrimination data.
%
% 


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

V1Mode = 6;
csfSelector = 'mean';

if ~exist('pars','var') || isempty(pars)
    pars = getPars(3,length(timecourse),V1Mode);
end

noiseConst = pars(1);
noiseFactor= pars(2);
pars(8:9) = round(pars(8:9));


if any(pars<0) 
    error('You passed negative parameters to the model. This is not allowed!');
end

[imageBitmap,degSize] = foveate_image0(imageBitmap,degSize,timecourse,csfSelector);

imageBitmap = imresize(imageBitmap,sizePx);

[outN,~,out] = V1(imageBitmap,degSize,4,pars(3:14),[]);

imageNoise = noiseConst+noiseFactor.*abs(out);
imageNoise = imageNoise./prod(degSize)*numel(imageBitmap).*pars(8)./12.*pars(9)./8;

