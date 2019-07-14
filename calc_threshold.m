function  [threshold,pCorrect] = calc_threshold(imSize,degSize,timecourse,freq,contrast0,orientation,L,winSize,phase,fixation,window,type,pars,tolerance,pedestalFreq,pedestalOrientation,pedestalPhase,sizePx,PC,V1Mode)
% function  calc_threshold(imSize,degSize,freq,contrast0,orientation,L,winSize,phase,fixation,window,type,pars,tolerance,pedestalFreq,pedestalOrientation,pedestalPhase,sizePx)
% this tries different contrasts to find a value with 75 % correct
% +-tolerance

if ~exist('tolerance','var') || isempty(tolerance)
    tolerance = 0.05;
end

if ~exist('pedestalFreq','var') || isempty(pedestalFreq)
    pedestalFreq = freq;
end
if ~exist('pedestalOrientation','var') || isempty(pedestalOrientation)
    pedestalOrientation = orientation;
end
if ~exist('pedestalPhase','var') || isempty(pedestalPhase)
    pedestalPhase = phase;
end
if ~exist('sizePx','var')
    sizePx = [];
end
if ~exist('PC','var') || isempty(PC)
    PC = .75;
end
if ~exist('V1Mode','var') || isempty(V1Mode)
    V1Mode = 7;
end
pTest= NaN;

if false %contrast0 == 0
    f0   = 4.1726;
    f1   = 1.3625;
    a    = 0.8493;
    p    = 0.7786;
    gain = 373.08;
    csf_handle = @(f) gain.*(sech((f./f0).^p)-a.*sech(f./f1));
    
    contrast = 1./csf_handle(freq)+contrast0;
    [diff, noise, p] = grating_experiment(imSize,degSize,timecourse,freq,contrast,contrast0,orientation,L,winSize,phase,fixation,window,type,pars,pedestalFreq,pedestalOrientation,pedestalPhase,sizePx,[],V1Mode);
    k = 0;
    while abs(p-0.75)>=tolerance
        if p >= 0.99
            contrast = contrast.*(.25.*exp(-k./50)+(1-exp(-k./50)));
        elseif p>PC
            contrast = contrast.*(.75.*exp(-k./50)+(1-exp(-k./50)));
        elseif p>0.55
            contrast = contrast.*(1.25.*exp(-k./50)+(1-exp(-k./50)));
        else
            contrast = contrast.*(2.*exp(-k./50)+(1-exp(-k./50)));
        end
        while contrast <=contrast0
            contrast = contrast+0.01;
        end
        k = k+1;
        [diff, noise, p] = grating_experiment(imSize,degSize,time,timecourse,freq,contrast,contrast0,orientation,L,winSize,phase,fixation,window,type,pars,pedestalFreq,pedestalOrientation,pedestalPhase,sizePx,[],V1Mode);
        %fprintf('%d: contrast = %d, pc = %d\n',k,contrast,p);
    end
    
    threshold = contrast;
    
else % for discrimination use bisection method
    c1 = 0;
    c2 = 1-contrast0;
    [diff, noise, p2] = grating_experiment(imSize,degSize,timecourse,freq,c2,contrast0,orientation,L,winSize,phase,fixation,window,type,pars,pedestalFreq,pedestalOrientation,pedestalPhase,sizePx,[],V1Mode);
    
    if p2<.75
        threshold = c2;
    else
        finished = 0;
        k = 0;
        while ~finished
            k = k+1;
            cTest = (c2+c1)/2;
            [diff, noise, pTest] = grating_experiment(imSize,degSize,timecourse,freq,cTest,contrast0,orientation,L,winSize,phase,fixation,window,type,pars,pedestalFreq,pedestalOrientation,pedestalPhase,sizePx,[],V1Mode);
            if pTest >PC
                c2 = cTest;
            elseif pTest <PC
                c1 = cTest;
            end
            if (c2-c1)/c1 < tolerance
                finished = 1;
                threshold = (c2+c1)/2;
            end 
            %fprintf('%d: contrast = %d, pc = %d\n',k,cTest,pTest);
        end
    end
    
    pCorrect = pTest;
end