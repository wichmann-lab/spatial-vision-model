function [diff, noise, p, pGrad] = detection_experiment(signal,background,degSize,timecourse,contrast,L,fixation,type,pars,sizePx,lumBorder,Gradidx,V1Mode)
%function [diff, noise, p] = detection_experiment(signal,background,degSize,time,timecourse,contrast,L,fixation,type,pars,sizePx,V1Mode)
% This function expects at least a background a signal and their retinal
% size in degrees. It then computes how well these can be distinguished.
% The simulation assumes that signal and background have the same
% timecourse. 
% The background is always processed as is by the model, i.e. should
% ideally be given in cd/m^2 units
% The signal is expected as a [-1 1] scaled contrast image. It is added to
% the background at [contrast]. If you want you can change the reference
% luminance by changing L. 


if ~exist('L','var') || isempty(L)
    L = mean(background(:));
end

if ~exist('pars','var') || isempty(pars)
    pars = [];
end

if ~exist('fixation','var') || isempty(fixation)
    fixation = degSize/2;
end

if ~exist('timecourse','var') || isempty(timecourse)
    timecourse = []; % set in model now!
end

if ~exist('sizePx','var') 
    sizePx = []; % set in model 
end
if ~exist('lumBorder','var') || isempty(lumBorder)
    lumBorder = [0,inf];
end
if ~exist('Gradidx','var') || isempty(Gradidx)
    Gradidx = []; 
end
if ~exist('V1Mode','var') || isempty(V1Mode)
    V1Mode = [];
end

diff = zeros(size(contrast));
noise = zeros(size(contrast));
p = zeros(size(contrast));

for icontrast = 1:length(contrast)
    blank = background;
    target = background + signal.*contrast.*L;
    target(target<lumBorder(1)) = lumBorder(1);
    target(target>lumBorder(2)) = lumBorder(2);
    blank(blank<lumBorder(1)) = lumBorder(1);
    blank(blank>lumBorder(2)) = lumBorder(2);
    if nargout >3
        [target,targetNoise,targetGrad,targetnoiseGrad] = early_vision_model(target,degSize,timecourse,fixation,pars,sizePx,Gradidx,[],V1Mode);
        [blank,blankNoise,blankGrad,blanknoiseGrad]   = early_vision_model(blank ,degSize,timecourse,fixation,pars,sizePx,Gradidx,[],V1Mode);
    else
        [target,targetNoise] = early_vision_model(target,degSize,timecourse,fixation,pars,sizePx,[],[],V1Mode);
        [blank,blankNoise]   = early_vision_model(blank ,degSize,timecourse,fixation,pars,sizePx,[],[],V1Mode);
    end
    if nargout>3
        [diff(icontrast), noise(icontrast), p(icontrast),pGrad] = compare_images(target,targetNoise,blank,blankNoise,1,targetGrad,targetnoiseGrad,blankGrad,blanknoiseGrad);
    else
        [diff(icontrast), noise(icontrast), p(icontrast)] = compare_images(target,targetNoise,blank,blankNoise,1);
    end
end
