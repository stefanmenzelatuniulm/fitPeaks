close all;
clear all;
clc;

%-------------SETTINGS-------------

TR = 30000; %ms
gamma13C = 10.7084*10^6; %gyromagnetic ration of 13C (Hz/T)
B0 = 3; %B0 of the scanner (T)
rf_center_ppm = 180.43; %Center of the RF pulse (ppm) %IST SCHON RELATIV ZU TMS
bw = 500; %Bandwidth of the RF pulse (Hz)
chemicalSpecies = "Lactate"; %Name(s) of the chemical species in the spectrum
scaleHzFit = 10; %Scaling factor for the x axis only in the context of the Lorentzian fit (CARE: very sensitive)
scalePpmFit = 0.085;
increaseSliderStepResolutionFactor = 32; %Increases the default slider resolution by this factor, relevant for phase correction
increasePhi1LimitsFactor = 5;
annotationXOffset = 0; %Offset of fit parameter annotation in X direction, if there is significant overlap with the plot
floorOffset = 0.015;
path = "C:\Users\menze\Desktop\Matlab\MR_Data\105"; %Path to data
attemptLorentzianFit = false; %Attempt to fit single Lorentzian
attemptFIDFit = false;

%-------------END OF SETTINGS-------------

%Reco-PC returns complex signal (not spectrum)
load(path+"\data.mat");

%Average over the number of averages 
Data = double(mean(squeeze(Data),2));

fig = figure('WindowState', 'maximized');
X_FID = linspace(1,TR,length(Data));
X_FID_t = transpose(linspace(1,TR,length(Data)));
plot(X_FID, abs(Data)/max(abs(Data)));
legend("FID", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");

if attemptFIDFit
    fitfunction="C+M0*exp(-x/T2s)";
    coeffs=["C" "M0" "T2s"];
    options=fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],'Upper',[1 1 inf],'StartPoint',[0 0 1]);
    fttype = fittype(fitfunction,coefficients=coeffs);
    ft=fit(X_FID_t,abs(Data)/max(abs(Data)),fttype,options);
    coeffvals = coeffvalues(ft);
    ci = confint(ft,0.95);
    str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',"T_2^*",coeffvals(3),ci(:,3));
    annotation('textbox',[0.53+annotationXOffset 0.69 0.2 0.2],'String',['Relevant fit coefficient with 95% confidence bounds: ', strtrim(str1+"   (ms)")],'EdgeColor','none',"FitBoxToText","on", "Color","r","FontSize",8);
    hold on;
    ax = gca;
    plot(ax,ft,"r");
    legend("FID", "Fit with $$C+M_0 e^{-\frac{t}{T_2^*}}$$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
end

title("FID of "+chemicalSpecies);
xlabel("$t$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
xlim([min(X_FID) max(X_FID)]);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
saveas(fig, path+chemicalSpecies+"_FID.fig");
saveas(fig, path+chemicalSpecies+"_FID.svg");
saveas(fig, path+chemicalSpecies+"_FID.png");

close all;

%Determine spectrum
spectrum = fftshift(fft(Data));

%Determine x Axis
f_ref = gamma13C*B0; %Hz
rf_center = (10^(-6))*rf_center_ppm*f_ref; %Frequency of RF pulse center (Hz) relative zu TMS -> Tobi : SO RECHNET MAN PPM IN HZ UM, (f-fref/fref) nur einmal beim Scanner bereits kalibriert
X_Hz_rel = linspace(-bw/2, bw/2, length(Data)); %Chemical shift / Offresonance (Hz) relative to center of RF pulse
X_ppm_rel = (X_Hz_rel/f_ref)*10^6+rf_center_ppm;
X_Sample = linspace(1, length(Data), length(Data));

%Automatic 0th order phase correction with respect to highest peak
[highestPeak, maxIndex] = max(abs(spectrum));
phaseAngleRad = angle(spectrum(maxIndex));
Data = Data*exp(-1i*phaseAngleRad);
spectrum = fftshift(fft(Data));

%Plot phase angle of spectrum 
fig = figure('WindowState', 'maximized');
plot(X_Sample, angle(spectrum)+pi); 
legend("Phase angle", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("Phase angle of the spectrum of "+chemicalSpecies);
xlim([min(X_Sample) max(X_Sample)]);
xlabel("Sample", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Phase angle (rad)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, path+chemicalSpecies+"_phase.fig");
saveas(fig, path+chemicalSpecies+"_phase.svg");
saveas(fig, path+chemicalSpecies+"_phase.png");

%Plot modulus of spectrum
fig = figure('WindowState', 'maximized');
plot(X_ppm_rel, floorOffset+(abs(spectrum)-median(abs(spectrum)))/max(abs(spectrum)-median(abs(spectrum)))); %Rough baseline correction
legend("Spectrum", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("Spectrum of "+chemicalSpecies, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
xlim([min(X_ppm_rel) max(X_ppm_rel)]);
xlabel("Chemical Shift (ppm)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, path+chemicalSpecies+"_ppm_Modulus.fig");
saveas(fig, path+chemicalSpecies+"_ppm_Modulus.svg");
saveas(fig, path+chemicalSpecies+"_ppm_Modulus.png");

close all;

%Determine 0th and 1st order phase correction
disp("Adjust 0th and 1st order phase correction");
S.fh = figure('WindowState', 'maximized');
S.ax = axes('unit', 'normalized', 'position', [0.05 0.15 0.9 0.8]);
S.X_ppm_rel = X_ppm_rel;
S.Data = Data;
S.phi0 = 0;
S.phi1 = 0;
sReal = real(spectrum);
S.p1 = plot(S.ax, X_ppm_rel, sReal);
hold on;
S.pivot = X_ppm_rel(maxIndex);
hold on;
Y = median(sReal)*ones(length(X_ppm_rel));
Y = [Y(1) Y(end)];
X = [X_ppm_rel(1) X_ppm_rel(end)];
S.p2 = plot(S.ax, X, Y, "Color", "m");
S.l = xline(S.pivot);
S.phi0Slider = uicontrol('style', 'slider', 'unit','normalized', 'position', [0.1 0.073 0.8 0.025], 'min', -pi, 'max', pi, 'value', 0, 'sliderstep',[0.01 0.1]/increaseSliderStepResolutionFactor, 'callback', {@SliderCB, 'phi0'}); 
txtphi0 = uicontrol('Style', 'text', 'unit', 'normalized', 'position', [0 0.073 0.1 0.025], 'String', '0th order phase correction', 'BackgroundColor', "White", 'HorizontalAlignment', 'Center');
S.phi1Slider = uicontrol('style','slide', 'unit', 'normalized', 'position', [0.1 0.048 0.8 0.025], 'min', 0, 'max', length(X_ppm_rel), 'value', 0, 'sliderstep',[0.01 0.1]/(increasePhi1LimitsFactor*increaseSliderStepResolutionFactor), 'callback', {@SliderCB, 'phi1'});
txtphi1 = uicontrol('Style','text', 'unit', 'normalized', 'position', [0 0.048 0.1 0.025], 'String', '1st order phase correction', 'BackgroundColor', "White", 'HorizontalAlignment', 'Center');  
S.pivotSlider = uicontrol('style','slide', 'unit', 'normalized', 'position', [0.1 0.023 0.8 0.025], 'min', min(X_ppm_rel), 'max', max(X_ppm_rel), 'value', X_ppm_rel(maxIndex), 'sliderstep',[0.01 0.1]/increaseSliderStepResolutionFactor, 'callback', {@SliderCB, 'pivot'});
txtPivot = uicontrol('Style','text', 'unit', 'normalized', 'position', [0 0.023 0.1 0.025], 'String', 'Pivot', 'BackgroundColor', "White", 'HorizontalAlignment', 'Center');   
update(S);
guidata(S.fh, S);
title("Click on sliders, then select pivot and adjust sliders with arrow keys for 0th and 1st order phase correction, then minimize figure and press Enter in Command Window");
legend("Real part of the spectrum", "Median of the real part of the spectrum $$\approx$$ noise floor", "Pivot", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest")
xlim([min(X_ppm_rel) max(X_ppm_rel)]);
xlabel("Chemical Shift (ppm)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14); 
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

pause;

phi0 = get(S.phi0Slider, "Value");
phi1 = round(get(S.phi1Slider, "Value"));
p = get(S.pivotSlider, "Value");

%Apply 0th and 1st order phase correction
[~, pivotIndex] = min(abs(S.X_ppm_rel-p));
tempData = circshift(S.Data, -phi1);
realSpectrumTemp = real(fftshift(fft(tempData)));
C = exp(-1i*(angle(realSpectrumTemp(pivotIndex))+phi0));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
sReal = real(fftshift(C.*fft(tempData)));

%Baseline correction and normalization -> Sinnvoll? Tobi fragen
sReal = sReal-median(sReal); %Wenn das noise level konsequent über 0 wäre, wäre doch der median entsprechend über 0? nicht mean verwenden, da peaks nicht mit einbezogen werden sollen
sReal = sReal/max(sReal);

%Plot and fit Lorentzian
fig = figure('WindowState', 'maximized');
ax = gca;
plot(ax, X_Hz_rel, sReal);
title("Spectrum of "+chemicalSpecies, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
fitfunction="(1/pi)*((FWHM/2)/((x/"+num2str(scaleHzFit)+"-x0/"+num2str(scaleHzFit)+")^2+(FWHM/2)^2))";
coeffs=["FWHM" "x0"];
options=fitoptions('Method','NonlinearLeastSquares','Lower',[-inf -inf],'Upper',[inf inf],'StartPoint',[10 mean(X_Hz_rel)]);
fttype = fittype(fitfunction,coefficients=coeffs);
X_Hz_t = transpose(X_Hz_rel);
if attemptLorentzianFit
    ft=fit(X_Hz_t,sReal,fttype,options);
    coeffs = coeffnames(ft);
    coeffvals = coeffvalues(ft);
    ci = confint(ft,0.95);
    str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{1},coeffvals(1)*scaleHzFit,ci(:,1)*scaleHzFit); %Skalierung passt so -> Lorentzian kann entsprechend umgestellt werden
    str2 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{2},coeffvals(2),ci(:,2)); %Skalierung passt so -> Lorentzian kann entsprechend umgestellt werden
    str3 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',"Therefore T_2^*",1000/(coeffvals(1)*scaleHzFit*pi),1000*ci(:,1)*scaleHzFit/(pi*(coeffvals(1)*scaleHzFit)^2));
    annotation('textbox',[0.53+annotationXOffset 0.69 0.2 0.2],'String',['Fit coefficients with 95% confidence bounds: ', strtrim(str1+"   (Hz)"), strtrim(strrep(str2,"x0","\omega_0")+"   (Hz)"), strtrim(str3+"   (ms)")],'EdgeColor','none',"FitBoxToText","on", "Color","r","FontSize",8);
    hold on;
    plot(ax,ft,"r");
    legend("Real part of the spectrum", "Fit with Lorentzian", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
else
    legend("Real part of the spectrum", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
end
xlim([min(X_Hz_rel) max(X_Hz_rel)]);
xlabel("Chemical Shift (Hz)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14); 
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, path+chemicalSpecies+"_Hz.fig");
saveas(fig, path+chemicalSpecies+"_Hz.svg");
saveas(fig, path+chemicalSpecies+"_Hz.png");

close all;

%Plot in ppm
fig = figure('WindowState', 'maximized');
ax = gca;
plot(ax, X_ppm_rel, sReal);
title("Spectrum of "+chemicalSpecies, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 16);
fitfunction="(1/pi)*((FWHM/2)/((x/"+num2str(scalePpmFit)+"-x0/"+num2str(scalePpmFit)+")^2+(FWHM/2)^2))";
coeffs=["FWHM" "x0"];
options=fitoptions('Method','NonlinearLeastSquares','Lower',[-inf -inf],'Upper',[inf inf],'StartPoint',[1 rf_center_ppm]);
fttype = fittype(fitfunction,coefficients=coeffs);
X_ppm_t = transpose(X_ppm_rel);
if attemptLorentzianFit
    ft=fit(X_ppm_t,sReal,fttype,options);
    coeffs = coeffnames(ft);
    coeffvals = coeffvalues(ft);
    ci = confint(ft,0.95);
    ci(:,1) = ci(:,1)*(10^(-6))*f_ref;
    coeffvals(1) = coeffvals(1)*(10^(-6))*f_ref;
    str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{1},coeffvals(1)*scalePpmFit,ci(:,1)*scalePpmFit); %Skalierung passt so -> Lorentzian kann entsprechend umgestellt werden
    str2 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{2},coeffvals(2),ci(:,2)); %Skalierung passt so -> Lorentzian kann entsprechend umgestellt werden
    str3 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',"Therefore T_2^*",1000/(coeffvals(1)*scalePpmFit*pi),1000*ci(:,1)*scalePpmFit/(pi*(coeffvals(1)*scalePpmFit)^2));
    annotation('textbox',[0.53+annotationXOffset 0.69 0.2 0.2],'String',['Fit coefficients with 95% confidence bounds: ', strtrim(str1+"   (Hz)"), strtrim(strrep(str2,"x0","\omega_0")+"   (ppm)"), strtrim(str3+"   (ms)")],'EdgeColor','none',"FitBoxToText","on", "Color","r","FontSize",8);
    hold on;
    plot(ax,ft,"r");
    legend("Real part of the spectrum", "Fit with Lorentzian", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
else
    legend("Real part of the spectrum", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
end
xlim([min(X_ppm_rel) max(X_ppm_rel)]);
xlabel("Chemical Shift (ppm)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14); 
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, path+chemicalSpecies+"_ppm.fig");
saveas(fig, path+chemicalSpecies+"_ppm.svg");
saveas(fig, path+chemicalSpecies+"_ppm.png");

close all;

%Callback functions for sliders
function SliderCB(aSlider, EventData, Param)
    S = guidata(aSlider); 
    S.(Param) = get(aSlider, 'Value');  
    guidata(aSlider, S); 
    update(S); 
end

function update(S)
    p = S.pivot;
    [~, pivotIndex] = min(abs(S.X_ppm_rel-p));
    phi0 = S.phi0;
    phi1 = round(S.phi1);
    tempData = circshift(S.Data, -phi1);
    realSpectrumTemp = real(fftshift(fft(tempData)));
    C = exp(-1i*(angle(realSpectrumTemp(pivotIndex))+phi0));                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    realSpectrum = real(fftshift(C.*fft(tempData)));
    set(S.p1, 'YData', realSpectrum); 
    set(S.l, "Value", S.pivot);
end