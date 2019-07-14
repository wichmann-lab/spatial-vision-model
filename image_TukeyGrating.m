function imageBitmap = image_TukeyGrating(imSize,degSize,freq,contrast,center,orientation,L,winSize,phase,alpha)
%function imageBitmap = image_hanningGrating(imSize,degSize,freq,contrast,center,orientation,L,winSize,phase,alpha)
% creates a sinwave grating with a Tukey window 

global useGPU
if isempty(useGPU)
    useGPU = false;
end
if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.5;
end

if ~exist('antiAliasing','var') || isempty(antiAliasing)
    antiAliasing = 4;
end
imSize = antiAliasing*imSize;


if useGPU
    x = gpuArray.linspace(0, degSize(2),imSize(2))-center(2);
    y = gpuArray.linspace(0, degSize(1),imSize(1))-center(1);
else
    x = linspace(0, degSize(2),imSize(2))-center(2);
    y = linspace(0, degSize(1),imSize(1))-center(1);
end

[xx,yy] = meshgrid(x,y');
xRot = sin(orientation)*yy+cos(orientation)*xx;
yRot = cos(orientation)*yy-sin(orientation)*xx;

if freq~=0
    grating = sin(2*pi*xRot*freq+phase);
else
    grating = 1+0*xRot;
end

xRot = xRot./winSize(1);
yRot = yRot./winSize(2);


window = image_window_tukey(imSize,degSize,center,winSize,orientation,alpha);

imageBitmap= L + contrast*L.*window.*grating;

imageBitmap = imresize(imageBitmap,1/antiAliasing);