function lumImage = get_lum_image(rgbImage,spectra)

if ~exist('spectra','var') || isempty(spectra)
    sfile = importdata('Color_Monitor.dat');
    spectra = sfile.data;
    clear sfile;
end

%% find power function for data

coeff = nan(3,4);
dataX = (0:255)';
%op = optimoptions('lsqnonlin');
%op.Display = 'off';
for iC = 1:4
    dataY = spectra(:,iC);
    %powerH = @(P) ((P(1).*dataX.^P(2)+min(min(spectra)))-dataY)./dataY;
    powerH = @(P) sum(((P(1).*dataX.^P(2)+min(min(spectra)))-dataY).^2);
    %powerH = @(P) (P(1).*dataX.^P(2)+min(min(spectra)))-dataY;
    %coeff(1:2,iC) = lsqnonlin(powerH,[0.01,2],[],[],op);
    coeff(1:2,iC) = fminsearch(powerH,[0.01,2]);
end
coeff(3,:) = min(min(spectra));

lumImage = double(rgbImage);
for iC = 1:3
    lumImage(:,:,iC) = coeff(1,iC).*lumImage(:,:,iC).^coeff(2,iC);
end

lumImage = sum(lumImage,3);
