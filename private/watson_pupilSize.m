function pupilSize = watson_pupilSize(luminance,age)
% function pupilSize = watson_pupilSize(image,cdMax.age)
% this function calculates the average light adapted pupil size according
% to Watson 2012
% This is necessary for the optics computations.
% Luminance should here be the "effective corneal flux density":
% L = L * area [deg.^2]
% times .1 if only one eye is used
% I included the more hypothetical formula for ages < 20

f = luminance.^.41;


if age >=20 && age <= 83
    pupilSize = (18.5172 + 0.122165.*f - 0.105569.*age +0.000138645.*f.*age) ...
        ./ (2+0.0630635.*f); % equation 27 (appendix 2)
elseif age >= 1 && age <= 83
    pupilSize = (16.4674 +exp(-0.208269.*age).*(-3.96868+0.00521209.*f)+0.124857.*f)...
        ./ (2+0.0630635.*f); % equation 28 (appendix 2)
else
    error('age out of well described bounds');
end

% Derivation for age smaller 20:
% Y = 19.291+44.03.*exp(-age./4.8441); %equivalent age
% pupilSize2 = (18.5172 + 0.122165.*f - 0.105569.*Y +0.000138645.*f.*Y) ...
%        ./ (2+0.0630635.*f); % equation 27 (appendix 2)
