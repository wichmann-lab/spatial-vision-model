function imageBitmap = image_hanningPlaid(imSize,degSize,freq,contrast,center,orientation,L,winSize,phase,antiAliasing)
%function imageBitmap = image_hanningGrating(imSize,degSize,freq,contrast,center,orientation,L,winSize,phase,antiAliasing)
% creates a sinwave plaid with the given parameters, you may pass a vector
% of two numbers in freq, contrast, orientation, and phase if you want to
% specify assymetric plaids, otherwise the values are mirrored

global useGPU
if isempty(useGPU)
    useGPU = false;
end

if ~exist('antiAliasing','var') || isempty(antiAliasing)
    antiAliasing = 4;
end
imSize = antiAliasing*imSize;

if useGPU
    imageBitmap = L.*gpuArray.ones(imSize);
else
    imageBitmap = L.*ones(imSize);
end

if ~all(contrast==0)
    if numel(freq) == 1
        freq = [freq,freq];
    end
    if numel(contrast) == 1
        contrast = [contrast,contrast];
    end
    if numel(orientation) == 1
        orientation = [pi/2+orientation,pi/2-orientation];
    end
    if numel(phase) == 1
        phase = [phase,phase];
    end
    assert(numel(freq)==numel(contrast) && numel(contrast) == numel(orientation) && numel(orientation) == numel(phase), 'Number of Contrasts, frequencies orientations and phases for the plaids should be equal or all in 1-2')
    
    if useGPU
        x = gpuArray.linspace(0, degSize(2),imSize(2))-center(2);
        y = gpuArray.linspace(0, degSize(1),imSize(1))-center(1);
    else
        x = linspace(0, degSize(2),imSize(2))-center(2);
        y = linspace(0, degSize(1),imSize(1))-center(1);
    end
    
    [xx,yy] = meshgrid(x,y');
    
    for iPlaid = 1:length(contrast)
        xRot = sin(orientation(iPlaid))*yy+cos(orientation(iPlaid))*xx;
        yRot = cos(orientation(iPlaid))*yy-sin(orientation(iPlaid))*xx;
        
        
        xRot1 = xRot./winSize(1);
        yRot1 = yRot./winSize(2);
        
        r2 = xRot1.^2+yRot1.^2;
        
        window = (r2<1).*(0.5+0.5*cos(pi*(sqrt(r2))));
        if freq(iPlaid)~=0
            grating = sin(2*pi*xRot*freq(iPlaid)+phase(iPlaid));
        else
            grating = 1+0*xRot;
        end
        imageBitmap= imageBitmap + contrast(iPlaid)*L.*window.*grating;
    end
end

imageBitmap = imresize(imageBitmap,1/antiAliasing);