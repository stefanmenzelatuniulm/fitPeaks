close all;
clear all;
clc;

%-------------SETTINGS-------------

TR = 60000; %ms
gamma13C = 10.7084*10^6; %gyromagnetic ration of 13C (Hz/T)
B0 = 3; %B0 of the scanner (T)
rf_center_ppm = 169.7; %Center of the RF pulse (ppm) %IST SCHON RELATIV ZU TMS
bw = 1200; %Bandwidth of the RF pulse (Hz)
chemicalSpecies = "Lactate"; %Name(s) of the chemical species in the spectrum
scaleHzFit = 13; %Scaling factor for the x axis only in the context of the Lorentzian fit (CARE: very sensitive)
scalePpmFit = 0.3;
increaseSliderStepResolutionFactor = 32; %Increases the default slider resolution by this factor, relevant for phase correction
increasePhi1LimitsFactor = 5;
annotationXOffset = 0; %Offset of fit parameter annotation in X direction, if there is significant overlap with the plot
path = "C:\Users\Stefan Menzel\Desktop\Matlab\MR_Data\13_03_2023\SpectrumLac\31"; %Path to data
attemptLorentzianFit = true; %Attempt to fit single Lorentzian
attemptFIDFit = true;

%-------------END OF SETTINGS-------------

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
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, path+chemicalSpecies+"_FID.fig");
saveas(fig, path+chemicalSpecies+"_FID.svg");

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
xlabel("Sample", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Phase angle (rad)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, path+chemicalSpecies+"_phase.fig");
saveas(fig, path+chemicalSpecies+"_phase.svg");

%Plot modulus of spectrum
fig = figure('WindowState', 'maximized');
plot(X_Hz_rel, (abs(spectrum)-median(abs(spectrum)))/max(abs(spectrum)-median(abs(spectrum)))); %Rough baseline correction
legend("Spectrum", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
title("Modulus of the spectrum of "+chemicalSpecies);
xlabel("Chemical Shift (Hz)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, path+chemicalSpecies+"_Hz_Modulus.fig");
saveas(fig, path+chemicalSpecies+"_Hz_Modulus.svg");

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
S.phi1Slider = uicontrol('style','slide', 'unit', 'normalized', 'position', [0.1 0.048 0.8 0.025], 'min', -increasePhi1LimitsFactor*pi, 'max', increasePhi1LimitsFactor*pi, 'value', 0, 'sliderstep',[0.01 0.1]/(increasePhi1LimitsFactor*increaseSliderStepResolutionFactor), 'callback', {@SliderCB, 'phi1'});
txtphi1 = uicontrol('Style','text', 'unit', 'normalized', 'position', [0 0.048 0.1 0.025], 'String', '1st order phase correction', 'BackgroundColor', "White", 'HorizontalAlignment', 'Center');  
S.pivotSlider = uicontrol('style','slide', 'unit', 'normalized', 'position', [0.1 0.023 0.8 0.025], 'min', min(X_ppm_rel), 'max', max(X_ppm_rel), 'value', X_ppm_rel(maxIndex), 'sliderstep',[0.01 0.1]/increaseSliderStepResolutionFactor, 'callback', {@SliderCB, 'pivot'});
txtPivot = uicontrol('Style','text', 'unit', 'normalized', 'position', [0 0.023 0.1 0.025], 'String', 'Pivot', 'BackgroundColor', "White", 'HorizontalAlignment', 'Center');   
update(S);
guidata(S.fh, S);
title("Click on sliders, then select pivot and adjust sliders with arrow keys for 0th and 1st order phase correction, then minimize figure and press Enter in Command Window");
legend("Real part of the spectrum", "Median of the real part of the spectrum $$\approx$$ noise floor", "Pivot", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest")
xlabel("Chemical Shift (ppm)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14); 
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

pause;

phi0 = get(S.phi0Slider, "Value");
phi1 = get(S.phi1Slider, "Value")/max(X_ppm_rel);
pivot = get(S.pivotSlider, "Value");

%Apply 0th and 1st order phase correction
X_ppm_rel = X_ppm_rel-pivot*phi1;
X_Hz_rel = X_Hz_rel-(10^(-6))*pivot*phi1*f_ref;
C = exp(-1i*(phi0-pivot*phi1));
sReal = real(fftshift(fft(C.*Data)));

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
    legend("Amplitude (a. u.)", "Fit with Lorentzian", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
else
    legend("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
end
xlabel("Chemical Shift (Hz)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14); 
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, path+chemicalSpecies+"_Hz.fig");
saveas(fig, path+chemicalSpecies+"_Hz.svg");

close all;

%Plot in ppm
fig = figure('WindowState', 'maximized');
ax = gca;
plot(ax, X_ppm_rel, sReal);
title("Spectrum of "+chemicalSpecies, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
fitfunction="(1/pi)*((FWHM/2)/((x/"+num2str(scalePpmFit)+"-x0/"+num2str(scalePpmFit)+")^2+(FWHM/2)^2))";
coeffs=["FWHM" "x0"];
options=fitoptions('Method','NonlinearLeastSquares','Lower',[-inf -inf],'Upper',[inf inf],'StartPoint',[1 round(mean(X_ppm_rel))]);
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
    legend("Amplitude (a. u.)", "Fit with Lorentzian", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
else
    legend("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
end
xlabel("Chemical Shift (ppm)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14); 
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, path+chemicalSpecies+"_ppm.fig");
saveas(fig, path+chemicalSpecies+"_ppm.svg");

close all;

%Callback functions for sliders
function SliderCB(aSlider, EventData, Param)
    S = guidata(aSlider); 
    S.(Param) = get(aSlider, 'Value');  
    guidata(aSlider, S); 
    update(S); 
end

function update(S)
    X = S.X_ppm_rel;
    p = S.pivot;
    phi1Temp = S.phi1/max(X);
    X = X-p*phi1Temp;
    phi0Temp = S.phi0;
    tempData = S.Data;
    C = exp(-1i*(phi0Temp-p*phi1Temp));
    realSpectrum = real(fftshift(fft(C.*tempData)));
    set(S.p1, 'YData', realSpectrum); 
    set(S.p1, 'XData',  X); 
    % pivotMin = get(S.pivotSlider, "Min");
    % pivotMax = get(S.pivotSlider, "Max");
    % if p>pivotMax
    %     set(S.pivotSlider, "Max", p)
    % end
    % if p<pivotMin
    %     set(S.pivotSlider, "Min", p)
    % end
    % set(S.pivotSlider, "Value", p);
    set(S.l, "Value", p);
end