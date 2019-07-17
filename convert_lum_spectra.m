function [imLum,imCones] = convert_lum_spectra(image,spectra,silent)
%function [imCones,imLum] = conesFromImage(image,spectra,silent)
% this function uses the measured spectra and the cone fundamentals from
% CVRL (Stockmans lab) to compute cone activations for each pixel in the
% image
% This version does this going through all 256 spectra-> requires a
% spectrum for each level. 
% some speed improvement implemented now!


if ~exist('spectra','var') || isempty(spectra)
    load('private/spectra256-2015-03-16.mat')
end
if ~exist('silent','var') || isempty(silent)
    silent = false;
end




image = double(image);
coneSens = importdata('private/ss10q_1.csv');
coneSens(isnan(coneSens))= -inf;
coneSens(:,2:4) = 10.^coneSens(:,2:4);


vlamb = importdata('private/logCIE2008v10q_1.csv');
vlambf = vlamb(:,1);
vlamb = 10.^(vlamb(:,2));



imCone1 = zeros(size(image,1),size(image,2),1,size(image,4));
imCone2 = zeros(size(image,1),size(image,2),1,size(image,4));
imCone3 = zeros(size(image,1),size(image,2),1,size(image,4));
imLum   = zeros(size(image,1),size(image,2),1,size(image,4));
for iRGB = 1:3
    freq = spectra{iRGB}(:,1);
    lumAbs = interp1(vlambf,vlamb,freq);
    lumAbs(isnan(lumAbs))=0;
    c1Sens =interp1(coneSens(:,1),coneSens(:,2),freq);
    c1Sens(isnan(c1Sens))=0;
    c2Sens =interp1(coneSens(:,1),coneSens(:,3),freq);
    c2Sens(isnan(c2Sens))=0;
    c3Sens =interp1(coneSens(:,1),coneSens(:,4),freq);
    c3Sens(isnan(c3Sens))=0;
    c1out = zeros(256,1);
    c2out = zeros(256,1);
    c3out = zeros(256,1);
    LumOut = zeros(256,1);
    for iValue = 0:255
        spectrum = spectra{iRGB}(:,iValue+2);
        c1out(iValue+1) = 4*sum(spectrum.*c1Sens);
        c2out(iValue+1) = 4*sum(spectrum.*c2Sens);
        c3out(iValue+1) = 4*sum(spectrum.*c3Sens);
        LumOut(iValue+1) = 4*sum(spectrum.*lumAbs);
    end
    imCone1 = imCone1 + c1out(image(:,:,iRGB,:)+1);
    imCone2 = imCone2 + c2out(image(:,:,iRGB,:)+1);
    imCone3 = imCone3 + c3out(image(:,:,iRGB,:)+1);
    imLum   = imLum + LumOut(image(:,:,iRGB,:)+1);
    if ~silent
        fprintf('finished channel %d\n',iRGB);
    end
end

imCones = cat(3,imCone1,imCone2,imCone3);

imLum = 683.002*imLum;