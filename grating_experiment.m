function [diff, noise, p, pGrad,noiseGrad,diffGrad] = grating_experiment(imSize,degSize,timecourse,freq,contrast,contrast0,orientation,L,winSize,phase,fixation,window,pars,pedestalFreq,pedestalOrientation,pedestalPhase,sizePx,Gradidx,V1Mode)
%function [diff, noise, p, pGrad,noiseGrad,diffGrad] = grating_experiment(imSize,degSize,timecourse,freq,contrast,contrast0,orientation,L,winSize,phase,fixation,window,pars,pedestalFreq,pedestalOrientation,pedestalPhase,sizePx,Gradidx,V1Mode)
%
% This simulates the two specified responses and computes their difference,
% weighted by the signal to noise ratio and returns the difference noise
% and probability correct.
%
% This will now also be able to calculate a gradient for the probability
% correct, if you ask for only one experiment and do not use other
% likelihoods.


if ~exist('orientation','var') || isempty(orientation)
    orientation = 0;
end

if ~exist('L','var') || isempty(L)
    L = 100;
end

if ~exist('winSize','var') || isempty(winSize)
    winSize = [2,2];
end

if ~exist('phase','var') || isempty(phase)
    phase = pi/2;
end

if ~exist('pars','var') || isempty(pars)
    pars = [];
end

if ~exist('fixation','var') || isempty(fixation)
    fixation = degSize/2;
end

if ~exist('window','var') || isempty(window)
    window = 'Gauss';
end

if ~exist('timecourse','var') || isempty(timecourse) % set in the model
    timecourse = []; 
end

if ~exist('pedestalFreq','var') || isempty(pedestalFreq)
    pedestalFreq = freq;
end
if ~exist('pedestalOrientation','var') || isempty(pedestalOrientation)
    pedestalOrientation = orientation;
end
if ~exist('pedestalPhase','var') || isempty(pedestalPhase) || isnan(pedestalPhase)
    pedestalPhase = phase;
end
if ~exist('sizePx','var') || isempty(sizePx)
    sizePx = [];
end
if ~exist('Gradidx','var') || isempty(Gradidx)
    Gradidx = [];
end
if ~exist('V1Mode','var') || isempty(V1Mode)
    V1Mode = [];
end

diff  = nan(length(freq),length(contrast));
noise = nan(length(freq),length(contrast));
p     = nan(length(freq),length(contrast));
for ifreq = 1:length(freq)
    for icontrast = 1:length(contrast)
        switch window % take care with winSize! different meaning for different windows
            case {'Gauss','gauss','Gabor','gabor'}
                gabor = image_gabor(imSize,degSize,freq(ifreq),contrast(icontrast),degSize/2,orientation,L,winSize(1),winSize(2),phase);
                blank = image_gabor(imSize,degSize,pedestalFreq(ifreq),contrast0,degSize/2,pedestalOrientation,L,winSize(1),winSize(2),pedestalPhase);
            case {'Hanning','2D Hanning','hanning','hanning 2D'}
                gabor = image_hanningGrating(imSize,degSize,freq(ifreq),contrast(icontrast),degSize/2,orientation,L,winSize,phase);
                blank = image_hanningGrating(imSize,degSize,pedestalFreq(ifreq),contrast0,degSize/2,pedestalOrientation,L,winSize,pedestalPhase);
            case {'Tukey.5','Tukey'}
                gabor = image_TukeyGrating(imSize,degSize,freq(ifreq),contrast(icontrast),degSize/2,orientation,L,winSize,phase,.5);
                blank = image_TukeyGrating(imSize,degSize,pedestalFreq(ifreq),contrast0,degSize/2,pedestalOrientation,L,winSize,pedestalPhase,.5);
            case {'logGabor','loggabor','log-Gabor','log-gabor'}
                gabor = image_hanningGrating(imSize,degSize,freq(ifreq),contrast(icontrast),degSize/2,orientation,L,winSize,phase);
                blank = image_hanningGrating(imSize,degSize,pedestalFreq(ifreq),contrast0,degSize/2,pedestalOrientation,L,winSize,pedestalPhase);
        end
        % NOW: Add gabor to blank interval!
        gabor = gabor+blank-L; 
        refLum = watson_filter_eye(ones(size(gabor)),degSize,4);
        refLum = L.*mean(refLum(:));
        
        if nargout >3
            [gabor,gaborNoise,gabGrad,gabNoiseGrad] = early_vision_model(gabor,degSize,timecourse,fixation,pars,sizePx,Gradidx,refLum,V1Mode);
            [blank,blankNoise,blankGrad,blankNoiseGrad] = early_vision_model(blank,degSize,timecourse,fixation,pars,sizePx,Gradidx,refLum,V1Mode);
        else
            [gabor,gaborNoise] = early_vision_model(gabor,degSize,timecourse,fixation,pars,sizePx,[],refLum,V1Mode);
            [blank,blankNoise] = early_vision_model(blank,degSize,timecourse,fixation,pars,sizePx,[],refLum,V1Mode);
        end
        if nargout > 3 % works only for one contrast and one frequency!
            [diff(ifreq,icontrast), noise(ifreq,icontrast), p(ifreq,icontrast),pGrad,noiseGrad,diffGrad] = compare_images(gabor,gaborNoise,blank,blankNoise,gabGrad,gabNoiseGrad,blankGrad,blankNoiseGrad);
        else
            [diff(ifreq,icontrast), noise(ifreq,icontrast), p(ifreq,icontrast)] = compare_images(gabor,gaborNoise,blank,blankNoise);
        end
    end
end
