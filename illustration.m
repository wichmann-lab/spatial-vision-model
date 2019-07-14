%% script with commands for illustration figures showing processing
% use the currently saved parameters
pars = getPars();

%% get original images
% A grating
grating = image_hanningGrating([512,512],[2,2],2,.1,[1,1],1,85,[1,1],1);
% A natural image
natImage = imread('DSC_0339.JPG');
natImageSmall=imresize(natImage,[96,96]);

%% show original image
image(natImage)
axis square
set(gca,'Visible','off')


%% convert to luminance and show again
[~,natImage] = conesFromImage(natImage);
natImageSmall=imresize(natImage,[96,96]);

natImageBackground = imresize(natImage,[256,256]);
natImageBackground = natImageBackground-min(natImageBackground(:));
natImageBackground = natImageBackground./max(natImageBackground(:));
natImageBackground = imadjust(natImageBackground,[min(natImageBackground(:)),max(natImageBackground(:))],[0,1],1/2.2);


% linear scale in luminance looks strange...
%imagesc(natImage)

% instead with gamma = 2.2
figure
J = natImage./max(natImage(:));
J = imadjust(J,[min(J(:)),max(J(:))],[0,1],1/2.2);
imagesc(J)
axis square
set(gca,'Visible','off')
colormap(gray(256))

figure
JG = grating./mean(grating(:))./2;
JG = imadjust(JG,[0,1],[0,1],1/2.2);
imagesc(JG,[0,1])
axis square
set(gca,'Visible','off')
colormap(gray(256))

%% Optics
natImage = watson_filter_eye(natImage,[2,2],4);
grating = watson_filter_eye(grating,[2,2],4);

figure
J = natImage./max(natImage(:));
J = imadjust(J,[min(J(:)),max(J(:))],[0,1],1/2.2);
imagesc(J)
axis square
set(gca,'Visible','off')
colormap(gray(256))

figure
JG = grating./mean(grating(:))./2;
JG = imadjust(JG,[0,1],[0,1],1/2.2);
imagesc(JG,[0,1])
axis square
set(gca,'Visible','off')
colormap(gray(256))
%% "fovea cut out

natImageCut = natImage./mean(natImage(:))-1;
gratingCut = grating./mean(grating(:))-1;

natImageCut = cutFovea(natImageCut,[2,2],[1.5,1.5],256);
gratingCut = cutFovea(gratingCut,[2,2],[1.5,1.5],256);

figure
J = (natImageCut-min(natImageCut(:)))./(max(natImageCut(:))-min(natImageCut(:)));
J = imadjust(J,[min(J(:)),max(J(:))],[0,1],1/2.2);
imagesc(J)
axis square
set(gca,'Visible','off')
colormap(gray(256))

figure
JG = (gratingCut-min(gratingCut(:)))./(max(gratingCut(:))-min(gratingCut(:)));
JG = imadjust(JG,[0,1],[0,1],1/2.2);
imagesc(JG,[0,1])
axis square
set(gca,'Visible','off')
colormap(gray(256))

%% full foveation

nframes = 228;
timecourse = 0.5 -0.5*cos((0:(nframes-1))/(nframes-1)*2*pi);
timecourse = imresize(timecourse,[1,1497],'nearest')';
[natImageLong,degSize]   = foveate_image1(natImage,[2,2],timecourse,[1,1],[],256);
[natImageMiddle,degSize] = foveate_image1(natImage,[2,2],ones(250,1),[1,1],[],256);
[natImageShort,degSize]  = foveate_image1(natImage,[2,2],ones(20,1),[1,1],[],256);

figure
for iPlot = 1:3
    subplot(1,3,iPlot)
    switch iPlot
        case 1
            imagesc(natImageLong,[-1.75,3])
        case 2
            imagesc(natImageMiddle,[-1.75,3])
        case 3
            imagesc(natImageShort,[-1.75,3])
    end
    axis square
    set(gca,'Visible','off')
    colormap(gray(256))
end

figure
imagesc(natImageMiddle)
axis square
set(gca,'Visible','off')
colormap(gray(1000))


%% V1 stuff
% parameters
noiseConst = pars(1);
noiseFactor = pars(2);
CNaka = pars(3);    % Naka rushton constant
ExNaka = pars(4);   % Naka rushton Exponent top
ExNakaNorm = pars(5); % Naka rushton Exponent below
bw(1) = pars(6);    % bandwidth in frequency (std of log-Gabor in octaves)
bw(2) = pars(7);    % bandwidth in orientation (std log-Gabor in radiants)
nFreq = pars(8);    % number of frequency bands
nOrient = pars(9);  % number of orientations
poolSize = pars(10); % meaning depends on type, spatial size of normalization pool
poolbw = pars(11);   % bandwidth of the normalization pool (Octaves in frequency)
minF   = pars(12);   % lowest frequency band
maxF   = pars(13);  % highest frequency band

%% Decomposition

[outNatImage,~,~,~,filters] = decomp_Gabor(natImageMiddle,degSize,[minF,maxF],nFreq,nOrient,bw);

figure
imagesc(flip(squeeze(mean(mean(abs(outNatImage),1),2)),2))
set(gca,'FontSize',14);
%xlabel('Frequency [cyc/deg]','FontSize',14)
%ylabel('Orientation [deg]','FontSize',14)
set(gca,'XTick',1:nFreq);
freq = 2.^linspace(log2(minF),log2(maxF),nFreq);
set(gca,'XTicklabel',sprintf('%.1f\n',freq))
orient = linspace(0,(nOrient-1)./nOrient*180,nOrient);
set(gca,'YTicklabel',sprintf('%.1f\n',orient))
%title('Raw Frequency Decomposition');
colormap(gray(256));
set(gca,'XTickLabelRotation',90)

figure
overlay = abs(outNatImage(:,:,5,7));
overlay = overlay./max(overlay(:));
im = cat(3,natImageBackground.*overlay,natImageBackground.*overlay,natImageBackground.*(1-overlay));
image(im)
axis square
set(gca,'Visible','off')

allNumbers = real(outNatImage(:,:,:,5));
for iChannel = 1:8
    natImageSmall=imresize(real(outNatImage(:,:,iChannel,5)),[96,96]);
    natImageSmall=natImageSmall-min(allNumbers(:));
%    imwrite(imadjust(natImageSmall./(max(allNumbers(:))-min(allNumbers(:))),[0,1],[0,1],1),sprintf('../../../figures/earlyVision/GifMaking/%03d.jpg',iChannel+15))
end


allNumbers = real(outNatImage(:,:,8,:));
for iChannel = 1:12
    natImageSmall=imresize(real(outNatImage(:,:,8,iChannel)),[96,96]);
    natImageSmall=natImageSmall-min(allNumbers(:));
%    imwrite(imadjust(natImageSmall./(max(allNumbers(:))-min(allNumbers(:))),[0,1],[0,1],1),sprintf('../../../figures/earlyVision/GifMaking/%03d.jpg',iChannel+15+8))
end
%% Normalization -> full V1

outNatImageNormalized = V1(natImageMiddle,degSize,1,pars(3:14));

figure
imagesc(flip(squeeze(mean(mean(outNatImageNormalized,1),2)),2))
set(gca,'FontSize',14);
%xlabel('Frequency [cyc/deg]','FontSize',14)
%ylabel('Orientation [deg]','FontSize',14)
set(gca,'XTick',1:nFreq);
freq = 2.^linspace(log2(minF),log2(maxF),nFreq);
set(gca,'XTicklabel',sprintf('%.1f\n',freq))
orient = linspace(0,(nOrient-1)./nOrient*180,nOrient);
set(gca,'YTicklabel',sprintf('%.1f\n',orient))
set(gca,'XTickLabelRotation',90)
%title('Activity after Normalization');
colormap(gray(256));

figure
imagesc(natImageBackground)
colormap(gray(256));
hold on
h =image(0.5*cat(3,ones(256),ones(256),ones(256)));
overlay = outNatImageNormalized(:,:,5,7);
set(h,'AlphaData',1-overlay./max(overlay(:)));
axis square
set(gca,'Visible','off')

figure
imagesc(natImageBackground)
colormap(gray(256));
hold on
h =image(cat(3,ones(256),ones(256),zeros(256)));
overlay = outNatImageNormalized(:,:,5,7);
set(h,'AlphaData',overlay./max(overlay(:)));
axis square
set(gca,'Visible','off')

figure
overlay = outNatImageNormalized(:,:,5,7);
overlay = overlay./max(overlay(:));
im = cat(3,natImageBackground.*overlay,natImageBackground.*overlay,natImageBackground.*(1-overlay));
image(im)
axis square
set(gca,'Visible','off')

%% show individual images

res =get(0,'ScreenPixelsPerInch');

h = figure;
pause(.5);
imagesc(real(outNatImage(:,:,5,12)))
colormap(gray(256));
hold on
set(gca,'Visible','off')
set(gca,'Units','pixels')
set(h,'Units','pixels')
set(h,'PaperUnits','inches')
set(h,'PaperPositionMode', 'manual');
set(gca,'Position',[1,1,257,257])
set(h,'PaperSize',[256,256]./res);
set(h,'PaperPosition',[0,0,256,256]./res)
pos = get(h,'Position');
set(h,'Position',[pos(1),pos(2),256,256])
drawnow
%print('../../../figures/EVIllustration/Channel1.jpg','-djpeg','-r96')

imagesc(real(outNatImage(:,:,4,8)))
colormap(gray(256));
drawnow
%print('../../../figures/EVIllustration/Channel2.jpg','-djpeg','-r96')

imagesc(real(outNatImage(:,:,1,5)))
colormap(gray(256));
drawnow
%print('../../../figures/EVIllustration/Channel3.jpg','-djpeg','-r96')


imagesc(real(outNatImageNormalized(:,:,5,12)))
colormap(gray(256));
drawnow
%print('../../../figures/EVIllustration/Channel1N.jpg','-djpeg','-r96')

imagesc(real(outNatImageNormalized(:,:,4,8)))
colormap(gray(256));
drawnow
%print('../../../figures/EVIllustration/Channel2N.jpg','-djpeg','-r96')

imagesc(real(outNatImageNormalized(:,:,1,5)))
colormap(gray(256));
drawnow
%print('../../../figures/EVIllustration/Channel3N.jpg','-djpeg','-r96')

hold off
imagesc(natImageBackground)
colormap(gray(256));
hold on
him =image(cat(3,ones(256),ones(256),zeros(256)));
overlay = abs(outNatImage(:,:,5,12));
set(him,'AlphaData',overlay./max(overlay(:)));
set(gca,'Visible','off')
drawnow
%print('../../../figures/EVIllustration/Overlay1.jpg','-djpeg','-r96')


imagesc(natImageBackground)
colormap(gray(256));
hold on
him =image(cat(3,ones(256),ones(256),zeros(256)));
overlay = abs(outNatImage(:,:,4,8));
set(him,'AlphaData',overlay./max(overlay(:)));
set(gca,'Visible','off')
drawnow
%print('../../../figures/EVIllustration/Overlay2.jpg','-djpeg','-r96')

imagesc(natImageBackground)
colormap(gray(256));
hold on
him =image(cat(3,ones(256),ones(256),zeros(256)));
overlay = abs(outNatImage(:,:,1,5));
set(him,'AlphaData',overlay./max(overlay(:)));
set(gca,'Visible','off')
drawnow
%print('../../../figures/EVIllustration/Overlay3.jpg','-djpeg','-r96')

imagesc(natImageBackground)
colormap(gray(256));
hold on
him =image(cat(3,ones(256),ones(256),zeros(256)));
overlay = outNatImageNormalized(:,:,5,12);
set(him,'AlphaData',overlay./max(overlay(:)));
set(gca,'Visible','off')
drawnow
%print('../../../figures/EVIllustration/Overlay1N.jpg','-djpeg','-r96')


imagesc(natImageBackground)
colormap(gray(256));
hold on
him =image(cat(3,ones(256),ones(256),zeros(256)));
overlay = outNatImageNormalized(:,:,4,8);
set(him,'AlphaData',overlay./max(overlay(:)));
set(gca,'Visible','off')
drawnow
%print('../../../figures/EVIllustration/Overlay2N.jpg','-djpeg','-r96')

imagesc(natImageBackground)
colormap(gray(256));
hold on
him =image(cat(3,ones(256),ones(256),zeros(256)));
overlay = outNatImageNormalized(:,:,1,5);
set(him,'AlphaData',overlay./max(overlay(:)));
set(gca,'Visible','off')
drawnow
%print('../../../figures/EVIllustration/Overlay3N.jpg','-djpeg','-r96')


%% Phase

h = figure;
pause(.5);
imagesc(real(outNatImage(:,:,4,5)))
colormap(gray(256));
hold on
set(gca,'Visible','off')
set(gca,'Units','pixels')
set(h,'Units','pixels')
set(h,'PaperUnits','inches')
set(h,'PaperPositionMode', 'manual');
set(gca,'Position',[1,1,257,257])
set(h,'PaperSize',[256,256]./res);
set(h,'PaperPosition',[0,0,256,256]./res)
pos = get(h,'Position');
set(h,'Position',[pos(1),pos(2),256,256])
drawnow

colormap(gray(256));
imagesc(real(outNatImage(:,:,4,5)))
drawnow
%print('../../../figures/EVIllustration/Real.jpg','-djpeg','-r96')


imagesc(imag(outNatImage(:,:,4,5)))
drawnow
%print('../../../figures/EVIllustration/Imag.jpg','-djpeg','-r96')


imagesc(abs(outNatImage(:,:,4,5)),[0,0.589])
drawnow
%print('../../../figures/EVIllustration/Absolute.jpg','-djpeg','-r96')

imagesc(natImageMiddle,[min(natImageMiddle(:)),max(natImageMiddle(:))])
drawnow
%print('../../../figures/EVIllustration/FinalNeuschwanstein.jpg','-djpeg','-r96')

hold off
Filt = fftshift(ifft2(filters{5}(:,:,4)));

imagesc(imag(Filt(65:192,65:192)),[-max(abs(Filt(:))),max(abs(Filt(:)))])
colormap(gray(256));
hold on
set(gca,'Visible','off')
set(gca,'Units','pixels')
set(gca,'Position',[1,1,257,257])
drawnow
%print('../../../figures/EVIllustration/FilterImag.jpg','-djpeg','-r96')

imagesc(real(Filt(65:192,65:192)),[-max(abs(Filt(:))),max(abs(Filt(:)))])
drawnow
%print('../../../figures/EVIllustration/FilterReal.jpg','-djpeg','-r96')
