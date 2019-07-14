function imageBitmap = image_hanningGrating(imSize,degSize,freq,contrast,center,orientation,L,winSize,phase,antiAliasing)
%function imageBitmap = image_hanningGrating(imSize,degSize,freq,contrast,center,orientation,L,winSize,phase,antiAliasing)
% creates a sinwave grating with a hanning window around it.

global useGPU
if isempty(useGPU)
    useGPU = false;
end

if ~exist('antiAliasing','var') || isempty(antiAliasing)
    antiAliasing = 4;
end
imSize = antiAliasing*imSize;

if contrast==0
    if useGPU
        imageBitmap = L.*gpuArray.ones(imSize);
    else
        imageBitmap = L.*ones(imSize);
    end
else
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
    
    r2 = xRot.^2+yRot.^2;
    
    window = (r2<1).*(0.5+0.5*cos(pi*(sqrt(r2))));
    
    imageBitmap= L + contrast*L.*window.*grating;
end

imageBitmap = imresize(imageBitmap,1/antiAliasing);