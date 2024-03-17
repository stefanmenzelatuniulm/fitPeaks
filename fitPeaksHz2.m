close all;
clear all;
clc;

%-------------SETTINGS-------------

TR = 60000; %ms
gamma13C = 10.7084*10^6; %gyromagnetic ration of 13C (Hz/T)
B0 = 3; %B0 of the scanner (T)
rf_center_ppm = 169.7; %Center of the RF pulse (ppm) 
bw = 1200; %Bandwidth of the RF pulse (Hz)
chemicalSpecies = "Alanine"; %Name(s) of the chemical species in the spectrum
scaleHzFit = 10; %Scaling factor for the x axis only in the context of the Lorentzian fit (CARE: very sensitive)
increaseSliderStepResolutionFactor = 32; %Increases the default slider resolution by this factor, relevant for phase correction
annotationXOffset = 0; %Offset of fit parameter annotation in X direction, if there is significant overlap with the plot
path = "C:\Users\menze\Desktop\Matlab\MR_Data\13_03_2023\SpectrumAla\33"; %Path to data
attemptLorentzianFit = true; %Attempt to fit single Lorentzian
attemptFIDFit = true;

%-------------END OF SETTINGS-------------

load(path+"\data.mat");

%Average over the number of averages 
Data = mean(squeeze(Data),2);

%Wieso ist der FID z.B. bei Lactat so komisch?! -> Tobi fragen, TR zu kurz -> Echo?
fig = figure('WindowState', 'maximized');
X_FID = linspace(1,TR,length(Data));
X_FID_t = transpose(linspace(1,TR,length(Data)));
plot(X_FID, abs(Data)/max(abs(Data)));
legend("FID");

if attemptFIDFit
    fitfunction="A+M0*exp(-x/T2s)";
    coeffs=["A" "M0" "T2s"];
    options=fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],'Upper',[1 1 inf],'StartPoint',[0 0 1]);
    fttype = fittype(fitfunction,coefficients=coeffs);
    ft=fit(X_FID_t,double(abs(Data)/max(abs(Data))),fttype,options);
    coeffvals = coeffvalues(ft);
    ci = confint(ft,0.95);
    str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',"T_2^*",coeffvals(3),ci(:,3));
    annotation('textbox',[0.53+annotationXOffset 0.69 0.2 0.2],'String',['Fit coefficients with 95% confidence bounds: ', strtrim(str1+"   (ms)")],'EdgeColor','none',"FitBoxToText","on", "Color","r","FontSize",8);
    hold on;
    ax = gca;
    plot(ax,ft,"r");
    legend("FID", "Fit with exponential decay");
end

title("FID of "+chemicalSpecies);
xlabel("$t$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, path+chemicalSpecies+"_FID.fig");
saveas(fig, path+chemicalSpecies+"_FID.svg");

close all;

%Determine spectrum
spectrum = fftshift(fft(Data));

%Determine x Axis
f_ref = gamma13C*B0; %Hz
rf_center = rf_center_ppm*f_ref+f_ref; %Absolute frequency of RF pulse center (Hz) -> Tobi fragen: passt das so
X_Hz_rel = linspace(-bw/2, bw/2, length(Data)); %Chemical shift (Hz) relative to RF pulse center

%Automatic 0th order phase correction with respect to highest peak
[highestPeak, maxIndex] = max(abs(spectrum));
phaseAngleRad = angle(spectrum(maxIndex));
Data = Data*exp(-1i*phaseAngleRad);
spectrum = fftshift(fft(Data));

%Determine 0th and 1st order phase correction
disp("Adjust 0th and 1st order phase correction");
S.fh = figure('WindowState', 'maximized');
S.ax = axes('unit', 'normalized', 'position', [0.05 0.15 0.9 0.8]);
S.angleFunction = @(x, phi0, phi1) phi0+phi1*x;
S.X = X_Hz_rel;
S.Data = Data;
S.phi0 = 0;
S.phi1 = 0;
S.p1 = plot(S.ax, X_Hz_rel, real(spectrum));
hold on;
Y = median(real(spectrum))*ones(length(X_Hz_rel));
Y = [Y(1) Y(end)];
X = [X_Hz_rel(1) X_Hz_rel(end)];
S.p2 = plot(S.ax, X, Y, "Color", "m");
update(S);
S.phi0Slider = uicontrol('style', 'slider', 'unit','normalized', 'position', [0.1 0.05 0.8 0.025], 'min', -pi, 'max', pi, 'value', 0, 'sliderstep',[0.01 0.1]/increaseSliderStepResolutionFactor, 'callback', {@SliderCB, 'phi0'}); 
txtphi0 = uicontrol('Style', 'text', 'unit', 'normalized', 'position', [0 0.05 0.1 0.025], 'String', '0th order phase correction', 'BackgroundColor', "White", 'HorizontalAlignment', 'Center');
S.phi1Slider = uicontrol('style','slide', 'unit', 'normalized', 'position', [0.1 0.025 0.8 0.025], 'min', -10*pi, 'max', 10*pi, 'value', 0, 'sliderstep',[0.01 0.1]/(10*increaseSliderStepResolutionFactor), 'callback', {@SliderCB, 'phi1'});
txtphi1 = uicontrol('Style','text', 'unit', 'normalized', 'position', [0 0.025 0.1 0.025], 'String', '1st order phase correction', 'BackgroundColor', "White", 'HorizontalAlignment', 'Center');   
guidata(S.fh, S);
title("Click on sliders, then adjust sliders with arrow keys for 0th and 1st order phase correction, then minimize figure and press Enter in Command Window");
legend("Real part of the spectrum", "Median of the real part of the spectrum ~ noise floor")
xlabel(string("Chemical Shift (Hz) relative to $\mathrm{RF_{center}}$ = "+num2str(rf_center/10^6)+" MHz or "+num2str(rf_center_ppm)+" ppm"), "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14); %Tobi fragen: Wert erscheint eine Größenordnung zu hoch
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

pause;

phi0 = get(S.phi0Slider, "Value");
phi1 = get(S.phi1Slider, "Value");

%Apply 0th and 1st order phase correction
Phi = transpose(phi0 + phi1*X_Hz_rel/max(X_Hz_rel));
realData = real(Data).*cos(Phi)-imag(Data).*sin(Phi); %https://pubmed.ncbi.nlm.nih.gov/28218950/
imagData = real(Data).*sin(Phi)+imag(Data).*cos(Phi);
Data = realData + 1i*imagData;
sReal = real(fftshift(fft(Data)));

%KANN SPEKTREN MIT MEHREREN PEAKS NICHT KORREKT PHASEN -> BASELINE
%CORRECTION FALSCH?
%Baseline correction and normalization -> Sinnvoll? Tobi fragen
sReal = sReal-median(sReal); %Wenn das noise level konsequent über 0 wäre, wäre doch der median entsprechend über 0? nicht mean verwenden, da peaks nicht mit einbezogen werden sollen
sReal = sReal/max(sReal);

%Plot and fit Lorentzian
fig = figure('WindowState', 'maximized');
ax = gca;
plot(ax, X_Hz_rel, sReal);
title("Spectrum of "+chemicalSpecies, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
fitfunction="(1/pi)*((FWHM/2)/((x/"+num2str(scaleHzFit)+"-x0/"+num2str(scaleHzFit)+")^2+(FWHM/2)^2))";
coeffs=["FWHM" "x0"];
options=fitoptions('Method','NonlinearLeastSquares','Lower',[-inf -inf],'Upper',[inf inf],'StartPoint',[10 round(mean(X_Hz_rel))]);
fttype = fittype(fitfunction,coefficients=coeffs);
X_Hz_t = transpose(X_Hz_rel);
if attemptLorentzianFit
    ft=fit(X_Hz_t,double(sReal),fttype,options);
    coeffs = coeffnames(ft);
    coeffvals = coeffvalues(ft);
    ci = confint(ft,0.95);
    str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{1},coeffvals(1)*scaleHzFit,ci(:,1)*scaleHzFit); %Skalierung passt so -> Lorentzian kann entsprechend umgestellt werden
    str2 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{2},coeffvals(2),ci(:,2)); %Skalierung passt so -> Lorentzian kann entsprechend umgestellt werden
    str3 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',"Therefore T_2^*",1000/(coeffvals(1)*scaleHzFit*pi),1000*ci(:,1)*scaleHzFit/(pi*(coeffvals(1)*scaleHzFit)^2));
    annotation('textbox',[0.53+annotationXOffset 0.69 0.2 0.2],'String',['Fit coefficients with 95% confidence bounds: ', strtrim(str1+"   (Hz)"), strtrim(strrep(str2,"x0","\omega_0")+"   (Hz)"), strtrim(str3+"   (ms)")],'EdgeColor','none',"FitBoxToText","on", "Color","r","FontSize",8);
    hold on;
    plot(ax,ft,"r");
    legend("Amplitude (a. u.)", "Fit with Lorentzian");
else
    legend("Amplitude (a. u.)");
end
xlabel(string("Chemical Shift (Hz) relative to $\mathrm{RF_{center}}$ = "+num2str(rf_center/10^6)+" MHz or "+num2str(rf_center_ppm)+" ppm"), "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14); %Tobi fragen: Wert erscheint eine Größenordnung zu hoch
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, path+chemicalSpecies+"_Hz.fig");
saveas(fig, path+chemicalSpecies+"_Hz.svg");

close all;

%Callback functions for sliders
function SliderCB(aSlider, EventData, Param)
    S = guidata(aSlider); 
    S.(Param) = get(aSlider, 'Value');  
    update(S); 
    guidata(aSlider, S); 
end

function update(S)
    angle = transpose(S.angleFunction(S.X/max(S.X), S.phi0, S.phi1));  
    realData = real(S.Data).*cos(angle)-imag(S.Data).*sin(angle); %https://pubmed.ncbi.nlm.nih.gov/28218950/
    imagData = real(S.Data).*sin(angle)+imag(S.Data).*cos(angle);
    Data = realData + 1i*imagData;
    realSpectrum = real(fftshift(fft(Data)));
    set(S.p1, 'YData', realSpectrum);  
end