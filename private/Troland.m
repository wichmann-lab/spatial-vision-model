function [imTroland] = Troland(imageBitmap,degSize,maxCD,pupilSize,age)
% function [imTroland] = cdToTroland(imageBitmap,maxCD,pupilSize)
% converts an image from display space into Trolands at the retina
% We assume the image is given as 0 = 0 and maxCD is the luminance which
% corresponds to 1 in the imageBitmap. Usually the maximal luminance output
% of the monitor.
% Then the image is transformed by pupil size into trolands.
% 
% We assume that pupil size depends only on the summed intensity of the
% image, e.g. does not depend on the spatial configuration
% You may provide an age of the subject to adjust the pupil size
% accordingly. If you don't we assume a 25 year old subject
% 
% this is now obsolete as I included this calculation into the filter eye
% function


if ~exist('age','var') || isempty(age)
    age = 25;
end
if ~exist('degSize','var') || isempty(degSize)
    degSize = getDegSize;
end

if ~exist('pupilSize','var') || isempty(pupilSize)
    % calculate "effective corneal flux density":
    lum       =  degSize(1).*degSize(2).*maxCD.*mean(imageBitmap(:)); 
    % use watson formula
    pupilSize = watson_pupilSize(lum,age);
    %pupilSize = 2; %mm % replace by Watson formula for asymptotic pupil-size
end



imTroland = imageBitmap.*maxCD.*pupilSize.^2.*pi./4;
