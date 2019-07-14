function pars = getPars(time,V1Mode)
%function pars = getPars(time,V1Mode)
% This returns the appropriate parameters for different models, observers
% and presentation times.
% time is assumed to map to the 4 presentation times used in the paper:
% 19,79,500,1497 ms
% V1Mode should be either 1 or 2 for local or mean normalization
% respectively.
% If no parameters were saved for the requested observer we default to the
% 'standard' observer for 79ms presentation time and local normalization.

if ~exist('time','var') || isempty(time)
    time = 79;
end
if ~exist('V1Mode','var') || isempty(V1Mode)
    V1Mode = 1;
end

switch time
    case {19,20}
        if V1Mode == 2
            load pars/pars1497_local.mat
            load pars/pars20_local.mat
            pars = [pars20,pars1497(6:end)];
        else
            load pars/pars20
            pars = pars20;
        end
    case 79
        if V1Mode == 2
            load pars/pars79_local.mat
            load pars/pars1497_local.mat
            pars = [pars79,pars1497(6:end)];
        else
            load pars/pars79
            pars = pars79;
        end
    case 500
        if V1Mode == 2
            load pars/pars1497_local.mat
            pars = pars1497;
            load pars/parsMF_local.mat
            pars(1:5) = parsMF;
        else
            load pars/parsMF
            pars = parsMF;
        end
    case 1497
        if V1Mode == 2
            load pars/pars1497_local
            pars = pars1497;
        else
            load pars/pars1497
            pars = pars1497;
        end
    otherwise
        if V1Mode == 2
            load pars/pars79_local.mat
            load pars/pars1497_local.mat
            pars = [pars79,pars1497(6:end)];
        else
            load pars/pars79
            pars = pars79;
        end
end
if ~exist('pars','var') || isempty(pars)
    pars = getPars();
end