function imageBitmap = image_grating(imSize,degSize,center,orientation,freq,phase)
%function imageBitmap = image_grating(imSize,degSize,center,orientation)
% creates a circular symmetric tukey window (raised cosine with a uniform
% area inbetween. 
% size gives the radius of the window (flat+cosinefalloff)
% alpha scales the proportion of the radius used for the cosine falloff
% center and size are both given in degrees


if ~exist('antiAliasing','var') || isempty(antiAliasing)
    antiAliasing = 4;
end
imSize = antiAliasing*imSize;

x = linspace(0, degSize(2),imSize(2))-center(2);
y = linspace(0, degSize(1),imSize(1))-center(1);

[xx,yy] = meshgrid(x,y');
xRot = sin(orientation)*yy+cos(orientation)*xx;
%yRot = cos(orientation)*yy-sin(orientation)*xx;

imageBitmap = sin(2*pi*xRot*freq+phase);

imageBitmap = imresize(imageBitmap,1/antiAliasing);