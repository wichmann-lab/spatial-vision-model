function MTF = watson_MTF(frequencies,pupil,age,pigment)
% function watson_MTF(frequencies,pupil,age,pigment)
% This calculates the MTF of the average corrected eye as described in
% Watson 2013
% you may omit the age and pigment inputs. Then they will be replaced by a
% worst case young eye (40 years, .16 pigmentation)

global useGPU
if isempty(useGPU)
    useGPU = false;
end

u0          = pupil.*pi.*10^6./555./180;
uhat        = frequencies./u0;
if useGPU
    MTF_diffrac = gpuArray.zeros(size(frequencies));
else
    MTF_diffrac = zeros(size(frequencies));
end
uhat2       = uhat(uhat<1);
MTF_diffrac(uhat<1) = 2./pi.*(acos(uhat2)-uhat2.*sqrt(1-uhat2.^2));
u1          = 21.95-5.512.*pupil+0.3922.*pupil.^2;
MTF         = (1+(frequencies./u1).^2).^(-0.62).*sqrt(MTF_diffrac);

% Scatter attenuation from IJspeert
if ~exist('age','var') || isempty(age)
    age     = 40;% still more or less no influence of age
end
if ~exist('pigment','var') || isempty(pigment)
    pigment = .16;% worst case: blue eye caucasian
end

MTF = (1-pigment)./(1+pigment.*(age./70).^4).*MTF;