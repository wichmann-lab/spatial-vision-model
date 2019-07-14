function [degSize,minValue] = getDegSize
% this sets a global degSize, thus all functions can get it from here if it
% is not supplied explicitly

degSize = [24.86058,31.07573];% values from range files in Potsdam

if nargout >1
    minValue = [0.8287,1.0359];
end