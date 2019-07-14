function imageBitmap = image_window_tukey(imSize,degSize,center,winSize,orientation,alpha)
%function imageBitmap = image_window_tukey(imSize,degSize,center,winSize,orientation,alpha)
% creates a tukey window (raised cosine with a uniform area inbetween). 
% winSize gives the radius of the window (flat+cosinefalloff)
% alpha scales the proportion of the radius used for the cosine falloff
% ([0,1])
% orientation gives the orientation of the windows first axis
% center and size are both winSize in degrees


if ~exist('antiAliasing','var') || isempty(antiAliasing)
    antiAliasing = 4;
end
imSize = antiAliasing*imSize;

x = linspace(0, degSize(2),imSize(2))-center(2);
y = linspace(0, degSize(1),imSize(1))-center(1);

[xx,yy] = meshgrid(x,y');
xRot = sin(orientation)*yy+cos(orientation)*xx;
yRot = cos(orientation)*yy-sin(orientation)*xx;

xRot = xRot./winSize(1);
yRot = yRot./winSize(2);

r2 = xRot.^2+yRot.^2;

r0 = (1-alpha).^2;
imageBitmap = (r2<r0)+(r2>r0).*(r2<1).*(0.5+0.5*cos(pi*(sqrt(r2)-(1-alpha))/alpha));


imageBitmap = imresize(imageBitmap,1/antiAliasing);