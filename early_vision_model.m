function [out,imageNoise,outGrad,noiseGrad] = early_vision_model(imageBitmap,degSize,timecourse,fixation,pars,sizePx,Gradidx,refLum,V1Mode,csfSelector)
%function [imageBitmap,imageNoise] = early_vision_model(imageBitmap,degSize,timecourse,fixation,pars,sizePx,Gradidx,refLum,V1Mode,csfSelector)
% This is the model in which only a few orientations enter the
% normalization pool. This should fit oblique masking additionally. 
%
% Using Gradidx you can specify which derivatives shall be computed. Gradidx is
% assumed to be sorted.



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
    pars = getPars(3,length(timecourse),V1Mode);
end
noiseConst = pars(1);
noiseFactor= pars(2);
pars(8:9) = round(pars(8:9));


if any(pars<0) 
    error('You passed negative parameters to the model. This is not allowed!');
end

imageBitmap = watson_filter_eye(imageBitmap,degSize,4);

[imageBitmap,degSize] = foveate_image1(imageBitmap,degSize,timecourse,fixation,refLum,sizePx,csfSelector);

if nargout <= 2
    out = V1(imageBitmap,degSize,V1Mode,pars(3:14));
else
    GradidxV1 = Gradidx(Gradidx>=3)-2;
    [out,outGrad] = V1(imageBitmap,degSize,V1Mode,pars(3:14),GradidxV1);
end

imageNoise = noiseConst+noiseFactor.*abs(out);
fac = (1./prod(degSize)*numel(imageBitmap).*pars(8)./12.*pars(9)./8);
imageNoise = imageNoise.*fac;
if nargout > 3
    if ismember(1,Gradidx) && ismember(2,Gradidx)
        noiseGrad = cat(5,ones(size(out)),out,outGrad.*noiseFactor).*fac;
    elseif ismember(1,Gradidx)
        noiseGrad = cat(5,ones(size(out)),outGrad.*noiseFactor).*fac;
    elseif ismember(2,Gradidx)
        noiseGrad = cat(5,out,outGrad.*noiseFactor).*fac;
    else
        noiseGrad = cat(5,outGrad.*noiseFactor).*fac;
    end
end
if nargout > 2
    if ismember(1,Gradidx) && ismember(2,Gradidx)
        outGrad = cat(5,zeros(size(out)),zeros(size(out)),outGrad);
    elseif ismember(1,Gradidx)
        outGrad = cat(5,zeros(size(out)),outGrad);
    elseif ismember(2,Gradidx)
        outGrad = cat(5,zeros(size(out)),outGrad);
    else
        %outGrad = outGrad;
    end
end
