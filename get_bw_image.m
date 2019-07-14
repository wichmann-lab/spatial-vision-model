function im_bw = get_bw_image(image,spectra)
%function im_bw = get_bw_image(image,spectra)
% this function uses the measured spectra and the cone fundamentals from
% CVRL (Stockmans lab) to compute cone activations for each pixel in the
% image
% It then finds the best approximation from the grayscale spectrum to
% recreate a black and white image.
% this first version does this going through all 256 spectra-> requires a
% spectrum for each level. 



if ~exist('spectra','var') || isempty(spectra)
    sfile = importdata('Color_Monitor.dat');
    spectra = sfile.data;
end
image = double(image);


vlamb = importdata('logCIE2008v10q_1.csv');
vlambf = vlamb(:,1);
vlamb = 10.^(vlamb(:,2));

imLum   = zeros(size(image,1),size(image,2));
for iRGB = 1:3
    for iValue = 0:255
        imLum(image(:,:,iRGB)==iValue)   = imLum(image(:,:,iRGB)==iValue) + spectra(iValue+1,iRGB);
        fprintf('finished %d of %d, channel %d\n',iValue,255,iRGB);
    end
end


%% build up black and white image
bw_values = spectra(:,4);

im_bw = zeros(size(imLum));
for iValue = 1:254
    im_bw((imLum>((bw_values(iValue+1-1)+bw_values(iValue+1))/2)) & (imLum<=((bw_values(iValue+1+1)+bw_values(iValue+1))/2))) = iValue;
    fprintf('Assembling Image: %d of %d\n',iValue,255);
end


