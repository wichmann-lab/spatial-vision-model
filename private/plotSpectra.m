function plotSpectra(spectra,spectra2)
% this function uses the measured spectra and the cone fundamentals from
% CVRL (Stockmans lab) to compute a luminance function to show to you in a
% figure
% This is meant as an overview to check whether the monitor measurements
% were useable

%coneSens = importdata('cone_fundamentals.csv');
%coneSens(isnan(coneSens))= -inf;
%coneSens(:,2:4) = 10.^coneSens(:,2:4);

%vlamb = importdata('v_lambda.dat');
%vlambf = vlamb(:,1);
%vlamb = 10.^(vlamb(:,2));

if iscell(spectra)
    gammaCurve   = zeros(256,4);
    for iRGB = 1:4
        freq = spectra{iRGB}(:,1);
        lumAbs = getVLambda(freq);
        lumAbs(isnan(lumAbs))=0;
        for iValue = 0:255
            spectrum = spectra{iRGB}(:,iValue+2);
            gammaCurve(iValue+1,iRGB) = 4*sum(spectrum.*lumAbs);
        end
    end
    gammaCurve = 683.002*gammaCurve*pi/2;
else
    gammaCurve = spectra;
end

if nargin>1
    if iscell(spectra2)
        gammaCurve2   = zeros(256,4);
        for iRGB = 1:4
            freq = spectra2{iRGB}(:,1);
            lumAbs = getVLambda(freq);
            lumAbs(isnan(lumAbs))=0;
            for iValue = 0:255
                spectrum = spectra2{iRGB}(:,iValue+2);
                gammaCurve2(iValue+1,iRGB) = 4*sum(spectrum.*lumAbs);
            end
        end
        gammaCurve2 = 683.002*gammaCurve2*pi/2;
    else
        gammaCurve2 = spectra2;
    end
end
%gammaCurve = 683.002*gammaCurve*pi; %times strange K_m constant and times pi for the steradiant thing...
figure;
subplot(2,3,1)
plot(0:255,gammaCurve(:,1),'r')
hold on
plot(0:255,gammaCurve(:,2),'g')
plot(0:255,gammaCurve(:,3),'b')
if nargin>1
    plot(0:255,gammaCurve2(:,1),'r--')
    plot(0:255,gammaCurve2(:,2),'g--')
    plot(0:255,gammaCurve2(:,3),'b--')
end
box off
set(gca,'TickDir','out');
xlabel('Input Value')
ylabel('luminance[cd/m.^2]')

subplot(2,3,2)
plot(0:255,gammaCurve(:,1)./max(gammaCurve(:,1)),'r')
hold on
plot(0:255,gammaCurve(:,2)./max(gammaCurve(:,2)),'g')
plot(0:255,gammaCurve(:,3)./max(gammaCurve(:,3)),'b')
if nargin>1
    plot(0:255,gammaCurve2(:,1)./max(gammaCurve2(:,1)),'r--')
    plot(0:255,gammaCurve2(:,2)./max(gammaCurve2(:,2)),'g--')
    plot(0:255,gammaCurve2(:,3)./max(gammaCurve2(:,3)),'b--')
end
box off
set(gca,'TickDir','out');
xlabel('Input Value')
ylabel('normalized Luminancen [0,1]')

subplot(2,3,3)
plot(0:255,gammaCurve(:,4),'k')
hold on
plot(0:255,sum(gammaCurve(:,1:3),2),'k--')
if nargin>1
    plot(0:255,gammaCurve2(:,4),'c')
    plot(0:255,sum(gammaCurve2(:,1:3),2),'c--')
end
box off
set(gca,'TickDir','out');
xlabel('Input Value')
ylabel('luminance[cd/m.^2]')


if iscell(spectra)
subplot(2,3,4)
plot(spectra{1}(:,1),spectra{1}(:,end),'r')
if nargin>1 && iscell(spectra2)
    hold on
    plot(spectra2{1}(:,1),spectra2{1}(:,end),'r--')
end
box off
set(gca,'TickDir','out');
xlabel('Wave Length [nm]')
ylabel('Spectral Radiance')
axis tight

subplot(2,3,5)
plot(spectra{2}(:,1),spectra{2}(:,end),'g')
if nargin>1 && iscell(spectra2)
    hold on
    plot(spectra2{2}(:,1),spectra2{2}(:,end),'g--')
end
box off
set(gca,'TickDir','out');
xlabel('Wave Length [nm]')
ylabel('Spectral Radiance')
axis tight

subplot(2,3,6)
plot(spectra{3}(:,1),spectra{3}(:,end),'b')
if nargin>1 && iscell(spectra2)
    hold on
    plot(spectra2{3}(:,1),spectra2{3}(:,end),'b--')
end
box off
set(gca,'TickDir','out');
xlabel('Wave Length [nm]')
ylabel('Spectral Radiance')
axis tight
end