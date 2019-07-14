function  threshold = calc_threshold_image(signal,background,degSize,timecourse,L,fixation,type,pars,PC,tolerance,lumBorder,V1Mode,csfSelector)
% function  threshold = calc_threshold_image(signal,background,degSize,timecourse,L,fixation,type,pars,PC,tolerance,lumBorder,V1Mode)
% this tries different contrasts +-tolerance to find a value with PC % correct
% 

if ~exist('fixation','var') || isempty(fixation)
    fixation = degSize/2;
end

if ~exist('PC','var') || isempty(PC)
    PC = 0.75;
end
if ~exist('tolerance','var') || isempty(tolerance)
    tolerance = 0.05;
end
if ~exist('lumBorder','var') 
    lumBorder = [];
end
if ~exist('V1Mode','var') 
    V1Mode = [];
end
if ~exist('csfSelector','var') 
    csfSelector = [];
end
c1 = 0;
c2 = 1;
%[diff, noise, p1] = grating_experiment(imSize,degSize,time,timcourse,freq,c1,contrast0,orientation,L,winSize,phase,fixation,window,type,pars);
%[diff, noise, p2] = grating_experiment(imSize,degSize,time,timcourse,freq,c2,contrast0,orientation,L,winSize,phase,fixation,window,type,pars);
[diff, noise, p2] = detection_experiment(signal,background,degSize,timecourse,c2,L,fixation,type,pars,[],lumBorder,[],V1Mode,csfSelector);
if p2<PC
    threshold = c2;
else
    finished = 0;
    k = 0;
    while ~finished
        k = k+1;
        cTest = (c2+c1)/2;
        [diff, noise, pTest] = detection_experiment(signal,background,degSize,timecourse,cTest,L,fixation,type,pars,[],lumBorder,[],V1Mode,csfSelector);
        if pTest >PC
            c2 = cTest;
        elseif pTest <PC
            c1 = cTest;
        end
        if (c2-c1)/c1 < tolerance
            finished = 1;
            threshold = (c2+c1)/2;
        end
    end
    %fprintf('%d: contrast = %d, pc = %d\n',k,cTest,pTest);
end

end