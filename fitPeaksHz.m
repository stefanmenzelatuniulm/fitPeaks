close all;
clear all;
clc;

%-------------SETTINGS-------------

gamma13C = 10.7084*10^6; %gyromagnetic ration of 13C (Hz/T)
B0 = 3; %B0 of the scanner (T)
rf_center_ppm = 169.7; %Center of the RF pulse (ppm) 
bw = 1200; %Bandwidth of the RF pulse (Hz)
chemicalSpecies = "Lactate"; %Name(s) of the chemical species in the spectrum
fitIncludeHeightLimit = -inf; %Points in the spectrum that are larger will be considered when fitting the Lorentzian
scaleHzFit = 10; %Scaling factor for the x axis only in the context of the Lorentzian fit (CARE: very sensitive)
increaseSliderStepResolutionFactor = 32; %Increases the default slider resolution by this factor, relevant for phase correction
annotationXOffset = 0; %Offset of fit parameter annotation in X direction, if there is significant overlap with the plot
path = "C:\Users\menze\Desktop\Matlab\MR_Data\13_03_2023\SpectrumLac\31"; %Path to data
attemptFit = true; %Attempt to fit single Lorentzian

%-------------END OF SETTINGS-------------

load(path+"\data.mat");

%Average over the number of averages and determine spectrum
Data = mean(squeeze(Data),2);
spectrum = fftshift(fft(Data));

f_ref = gamma13C*B0; %Hz
rf_center = rf_center_ppm*f_ref+f_ref; %Absolute frequency of RF pulse center (Hz) -> Tobi fragen: passt das so
X_Hz_rel = linspace(-bw/2, bw/2, length(Data)); %Chemical shift (Hz) relative to RF pulse center

%Determine 0th order phase correction
disp("Adjust 0th order phase correction");
figure('WindowState', 'maximized');
ax = gca;
plotHandle = plot(ax, X_Hz_rel, real(spectrum));
hold on;
Y = median(real(spectrum))*ones(length(X_Hz_rel));
Y = [Y(1) Y(end)];
X = [X_Hz_rel(1) X_Hz_rel(end)];
plot(ax, X, Y, "Color", "m");
title("Click on slider, then adjust slider with arrow keys for 0th order phase correction, then minimize figure and press Enter in Command Window");
legend("Real part of the spectrum", "Median of the real part of the spectrum ~ noise floor")
xlabel(string("Chemical Shift (Hz) relative to $\mathrm{RF_{center}}$ = "+num2str(rf_center/10^6)+" MHz"), "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14); %Tobi fragen: Wert erscheint eine Größenordnung zu hoch
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
h = uicontrol('style','slider','units','pixel','position',[20 20 300 20],'Min',0,'Max',2*pi); %-> Tobi fragen: Grenzen sinnvoll? Kann nicht mehr als phase 2*pi von imag=0 entfernt sein
sliderStep = get(h, "SliderStep");
set(h, "SliderStep", sliderStep/increaseSliderStepResolutionFactor);
%h(2) = uicontrol('style','slider','units','pixel','position',[20 40 300
%20],'Min',0,'Max',pi); %NUR DIE ERSTE KOMPONENTE WIRD AN makeplot ÜBERGEBEN !/(%(!/"%/("§%
addlistener(h,'ContinuousValueChange', @(hObject, event) makeplot0(hObject, event, Data, plotHandle));
pause;
phi0 = get(h, "Value");

%Apply 0th order phase correction -> Tobi fragen: korrekt?
Data = Data.*exp(-1i*phi0);

close all;

%Select peak where 0th order phase correction was applied as pivot
disp("Select peak where 0th order phase correction was applied");
figure('WindowState', 'maximized');
ax = gca;
spectrum = fftshift(fft(Data));
sReal = real(spectrum);
plot(ax, X_Hz_rel, sReal);
hold on;
[mn, index] = min(abs(X_Hz_rel-mean(X_Hz_rel)));
plotHandle = plot(ax, mn, sReal(index),"Marker","o","MarkerSize",12);
title("Click on slider, then adjust slider with arrow keys to select peak where 0th order phase correction was applied, then minimize figure and press Enter in Command Window");
legend("Real part of the spectrum");
xlabel(string("Chemical Shift (Hz) relative to $\mathrm{RF_{center}}$ = "+num2str(rf_center/10^6)+" MHz"), "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14); 
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
h = uicontrol('style','slider','units','pixel','position',[20 20 300 20],'Min',min(X_Hz_rel),'Max',max(X_Hz_rel)); 
sliderStep = get(h, "SliderStep");
set(h, "SliderStep", [1/length(X_Hz_rel) 1/length(X_Hz_rel)]);
addlistener(h,'ContinuousValueChange', @(hObject, event) makeplot3(hObject, event, sReal, X_Hz_rel, plotHandle));
pause;
X0 = get(h, "Value");

%Determine 1st order phase correction
disp("Adjust 1st order phase correction");
figure('WindowState', 'maximized');
ax = gca;
plotHandle = plot(ax, X_Hz_rel, sReal);
hold on;
Y = median(real(spectrum))*ones(length(X_Hz_rel));
Y = [Y(1) Y(end)];
X = [X_Hz_rel(1) X_Hz_rel(end)];
plot(ax, X, Y, "Color", "m");
title("Click on slider, then adjust slider with arrow keys for 1st order phase correction, then minimize figure and press Enter in Command Window");
legend("Real part of the spectrum", "Median of the real part of the spectrum ~ noise floor")
xlabel(string("Chemical Shift (Hz) relative to $\mathrm{RF_{center}}$ = "+num2str(rf_center/10^6)+" MHz"), "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14); %Tobi fragen: Wert erscheint eine Größenordnung zu hoch
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
h = uicontrol('style','slider','units','pixel','position',[20 20 300 20],'Min',0,'Max',4*pi); %Welche Grenzen hat phi1, wenn X normiert ist -> Tobi fragen
minorTicks = get(h, "SliderStep");
set(h, "SliderStep", sliderStep/(4*increaseSliderStepResolutionFactor));
addlistener(h,'ContinuousValueChange', @(hObject, event) makeplot1(hObject, event, Data, X0, plotHandle));
pause;
phi1 = get(h, "Value");

%Apply 1st order phase correction -> Tobi fragen: korrekt?
Data = Data.*exp(-1i*phi1*(X_Hz_rel-X0)/max(X_Hz_rel));

close all;

%Real part of corrected spectrum
sReal = real(fftshift(fft(Data)));

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
indices = find(sReal>fitIncludeHeightLimit);
if attemptFit
    ft=fit(X_Hz_t(indices),double(sReal(indices)),fttype,options);
    coeffs = coeffnames(ft);
    coeffvals = coeffvalues(ft);
    ci = confint(ft,0.95);
    str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{1},coeffvals(1)*scaleHzFit,ci(:,1)*scaleHzFit);
    str2 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{2},coeffvals(2),ci(:,2)); %dont scale with scaleHzFit for some reason -> Tobi fragen
    str3 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',"Therefore T_2^*",1000/(coeffvals(1)*scaleHzFit*pi),1000*ci(:,1)*scaleHzFit/(pi*(coeffvals(1)*scaleHzFit)^2));
    annotation('textbox',[0.53+annotationXOffset 0.69 0.2 0.2],'String',['Fit coefficients with 95% confidence bounds: ', strtrim(str1+"   (Hz)"), strtrim(strrep(str2,"x0","\omega_0")+"   (Hz)"), strtrim(str3+"   (ms)")],'EdgeColor','none',"FitBoxToText","on", "Color","r","FontSize",8);
    hold on;
    plot(ax,ft,"r");
    hold on;
    plot(ax, X_Hz_t(indices),double(sReal(indices)),"o","Color","m");
    legend("Amplitude (a. u.)", "Fit with Lorentzian", "Points included into fit");
else
    legend("Amplitude (a. u.)");
end
xlabel(string("Chemical Shift (Hz) relative to $\mathrm{RF_{center}}$ = "+num2str(rf_center/10^6)+" MHz"), "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14); %Tobi fragen: Wert erscheint eine Größenordnung zu hoch
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

saveas(fig, path+chemicalSpecies+"_Hz.fig");
saveas(fig, path+chemicalSpecies+"_Hz.svg");

close all;

%Callback method for 0th order phase correction
function makeplot0(h, event, Data, plotHandle)
    phi0 = get(h, 'Value');
    set(plotHandle, 'ydata', real(fftshift(fft(Data.*exp(-1i*phi0)))));
    drawnow;
end

%Callback method for 1st order phase correction
function makeplot1(h, event, Data, X0, plotHandle)
    phi1 = get(h, 'Value');
    X = transpose(get(plotHandle, 'Xdata'));
    maxX = max(X);
    X = X/maxX;
    X0 = X0/maxX;
    DataCorr = Data.*exp(-1i*phi1*(X-X0));
    set(plotHandle, 'Ydata', real(fftshift(fft(DataCorr))));
    drawnow;
end

%Callback method for pivot selection after 0th order phase correction
function makeplot3(h, event, sReal, X, plotHandle)
    C = get(h, 'Value');
    [mn, index] = min(abs(X-C));
    set(plotHandle, "Xdata", X(index));
    set(plotHandle, "Ydata", sReal(index));
    drawnow;
end