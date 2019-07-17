function image = image_gabor(imSize,degSize,freq,contrast,center,orientation,L,sigma0,sigma1,phase,antiAliasing)
% function image = image_gabor(imSize,degSize,freq,contrast,center,orientation,L,sigma0,sigma1,phase,antiAliasing)
% this function creates a Gabor patch with given orientation, frequency and
% size (rad and deg of visual angle) in the center of the image


if ~exist('antiAliasing','var') || isempty(antiAliasing)
    antiAliasing = 4;
end
imSize = antiAliasing*imSize;

x = ((1:imSize(2))./imSize(2).*degSize(2)-center(2));
y = ((1:imSize(1))./imSize(1).*degSize(1)-center(1));

[xx,yy] = meshgrid(x,y');
xRot = sin(orientation)*yy+cos(orientation)*xx;
yRot = cos(orientation)*yy-sin(orientation)*xx;

grating = sin(2*pi*xRot*freq+phase);

gauss = exp(-0.5*(xRot.^2./sigma0.^2+yRot.^2./sigma1.^2));

image = L + contrast.*grating.*gauss.*L;


image = imresize(image,1/antiAliasing);