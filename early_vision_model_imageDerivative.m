function [out,imageNoise,dImage] = early_vision_model3_imageDerivative(imageBitmap,degSize,timecourse,fixation,pars,sizePx,refLum,V1Mode,csfSelector)
%function [imageBitmap,imageNoise] = early_vision_model3(imageBitmap,degSize,timecourse,fixation,pars,sizePx,Gradidx,refLum,V1Mode)
% This is the model in which only a few orientations enter the
% normalization pool. This should fit oblique masking additionally. 
%
% According to Itti, Koch & Brown (2000) this should also fit orientation
% and frequency discrimination data.



if ~exist('degSize','var') || isempty(degSize)
    degSize = getDegSize;
end

if ~exist('fixation','var') || isempty(fixation)
    fixation = degSize/2;
end
if ~exist('sizePx','var') || isempty(sizePx)
    sizePx = 256;
end
if ~exist('timecourse','var') || isempty(timecourse)
    timecourse = 80;
end
if ~exist('refLum','var') || isempty(refLum)
    refLum = [];
end
if isscalar(timecourse)
    timecourse = ones(ceil(timecourse),1);  %ms
end
if ~exist('V1Mode','var') || isempty(V1Mode)
    V1Mode = 6;
end
if ~exist('csfSelector','var') || isempty(csfSelector)
    if V1Mode == 6
        csfSelector = 'mean';
    elseif V1Mode == 7
        csfSelector = 'standard';
    else
        csfSelector = [];
    end
end

if ~exist('pars','var') || isempty(pars)
    pars = getPars(3,length(timecourse),V1Mode);
end
noiseConst = pars(1);
noiseFactor= pars(2);
pars(8:9) = round(pars(8:9));


if any(pars<0) 
    error('You passed negative parameters to the model. This is not allowed!');
end

[imageBitmap,degSize] = foveate_image1(imageBitmap,degSize,timecourse,fixation,refLum,sizePx,csfSelector);

if nargout <= 2
    out = V1(imageBitmap,degSize,V1Mode,pars(3:14));
else
    out = V1(imageBitmap,degSize,V1Mode,pars(3:14));
end

imageNoise = noiseConst+noiseFactor.*abs(out);
fac = (1./prod(degSize)*numel(imageBitmap).*pars(8)./12.*pars(9)./8);
imageNoise = imageNoise.*fac;

% call V1 image Derivative here
