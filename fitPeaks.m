close all;
clear all;
clc;

%-------------SETTINGS-------------

gamma13C = 10.7084*10^6; %Hz/T
B0 = 3; %T
rf_center_ppm = 169.7; 
bw = 1200; 
chemicalSpecies = "Lactate";
%corrOffset = 1;
fitIncludeHeightLimit = 0.05;
scaleHzFit = 10;
scalePpmFit = 1;
path = "C:\Users\menze\Desktop\Matlab\MR_Data\13_03_2023\SpectrumLac\31";

%-------------END OF SETTINGS-------------

load(path+"\data.mat");

Data = mean(squeeze(Data),2);
Data = phaseCorrect(Data);
spectrum = fft(Data);

f_ref = gamma13C*B0; %Hz
rf_center = rf_center_ppm*f_ref+f_ref; %Hz %Tobi fragen: passt das so
X_Hz = linspace(-bw/2+rf_center, bw/2+rf_center, length(Data)); %absolute frequency (Hz)
X_Hz_rel = linspace(-bw/2, bw/2, length(Data)); %Chemical shift (Hz)
X_ppm = ((X_Hz-f_ref)/f_ref-rf_center_ppm)*10^6; %Chemical shift (ppm) %Tobi fragen: wenn relativ to TMS, dann -f_TMS_ppm statt -rf_center_ppm

% 0th order phase correction with respect to highest peak
% [highestPeak, maxIndex] = max(abs(spectrum));
% phaseAngleRad = angle(spectrum(maxIndex+corrOffset));
% spectrum = spectrum*exp(-1i*phaseAngleRad);

%Baseline correction and normalization
sReal = real(spectrum)-median(real(spectrum));
sReal = sReal/max(sReal);

%TODO: 1st order phase correction: plot real{spectrum*exp[-1i*(X-X_0thOrderCorrected)*C]}, adapt C with slider
%until spectrum is optimized.
%Not necessary if every peak is acquired individually and 0th order phase
%corrected individually!

fig1 = figure('WindowState', 'maximized');
ax1 = gca;
plot(ax1, sReal);
xlabel("Sample", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
title("Spectrum of "+chemicalSpecies, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
saveas(fig1, path+chemicalSpecies+"_Samples.fig");
saveas(fig1, path+chemicalSpecies+"_Samples.svg");

fig2 = figure('WindowState', 'maximized');
ax2 = gca;
plot(ax2, X_Hz_rel, sReal);
title("Spectrum of "+chemicalSpecies, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
fitfunction="(1/pi)*((FWHM/2)/((x/"+num2str(scaleHzFit)+"-x0/"+num2str(scaleHzFit)+")^2+(FWHM/2)^2))";
coeffs=["FWHM" "x0"];
options=fitoptions('Method','NonlinearLeastSquares','Lower',[-inf -inf],'Upper',[inf inf],'StartPoint',[10 round(mean(X_Hz_rel))]);
fttype = fittype(fitfunction,coefficients=coeffs);
X_Hz_t = transpose(X_Hz_rel);
indices = find(sReal>fitIncludeHeightLimit);
ft=fit(X_Hz_t(indices),double(sReal(indices)),fttype,options);
coeffs = coeffnames(ft);
coeffvals= coeffvalues(ft);
ci = confint(ft,0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{1},coeffvals(1)*scaleHzFit,ci(:,1)*scaleHzFit);
str2 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{2},coeffvals(2)*scaleHzFit,ci(:,2)*scaleHzFit);
annotationXOffset = 0;
annotation('textbox',[0.53+annotationXOffset 0.69 0.2 0.2],'String',['Fit coefficients with 95% confidence bounds: ', str1+" (Hz)", str2+" (Hz)"],'EdgeColor','none',"FitBoxToText","on", "Color","r","FontSize",8);
hold on;
plot(ax2,ft,"r");
hold on;
plot(ax2, X_Hz_t(indices),double(sReal(indices)),"o","Color","m");
legend("Amplitude (a. u.)", "Fit with Lorentzian", "Points included into fit");
xlabel(string("Chemical Shift (Hz) relative to RFcenter = "+num2str(rf_center/10^6)+" MHz"), "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
saveas(fig2, path+chemicalSpecies+"_Hz.fig");
saveas(fig2, path+chemicalSpecies+"_Hz.svg");

fig3 = figure('WindowState', 'maximized');
ax3 = gca;
plot(ax3,X_ppm,sReal);
fitfunction="(1/pi)*((FWHM/2)/((x/"+num2str(scalePpmFit)+"-x0/"+num2str(scalePpmFit)+")^2+(FWHM/2)^2))";
coeffs=["FWHM" "x0"];
options=fitoptions('Method','NonlinearLeastSquares','Lower',[-inf -inf],'Upper',[inf inf],'StartPoint',[1 round(mean(X_Hz_rel))]);
fttype = fittype(fitfunction,coefficients=coeffs);
X_ppm_t = transpose(X_ppm);
ft=fit(X_ppm_t(indices),double(sReal(indices)),fttype,options);
coeffs = coeffnames(ft);
coeffvals= coeffvalues(ft);
ci = confint(ft,0.95);
str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{1},coeffvals(1)*scalePpmFit,ci(:,1)*scalePpmFit);
str2 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{2},coeffvals(2)*scalePpmFit,ci(:,2)*scalePpmFit);
annotationXOffset = 0;
annotation('textbox',[0.53+annotationXOffset 0.69 0.2 0.2],'String',['Fit coefficients with 95% confidence bounds: ', str1+" (ppm)", str2+" (ppm)"],'EdgeColor','none',"FitBoxToText","on", "Color","r","FontSize",8);
hold on;
plot(ax3,ft,"r");
title("Spectrum of "+chemicalSpecies, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
hold on;
plot(ax3, X_ppm_t(indices),double(sReal(indices)),"o","Color","m");
legend("Amplitude (a. u.)", "Fit with Lorentzian", "Points included into fit");
xlabel(string("Chemical shift (ppm) relative to RFcenter = "+num2str(rf_center_ppm)+" ppm"), "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
saveas(fig3, path+chemicalSpecies+"_ppm.fig");
saveas(fig3, path+chemicalSpecies+"_ppm.svg");

fig4 = figure('WindowState', 'maximized');
ax4 = gca;
plot(ax4, abs(spectrum));
xlabel("Samples", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
ylabel("Amplitude (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
title("Modulus of the spectrum of "+chemicalSpecies, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
saveas(fig4, path+chemicalSpecies+"_abs.fig");
saveas(fig4, path+chemicalSpecies+"_abs.svg");

%close all; 