function [] = amp_spec_plot(pltdata,samprate)
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.
%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
%---------------------------------------
%
% average amp spectrum of some data
%
%--------------------------------------
%
nt = length(pltdata);
dt = samprate/1.e6;

% invariants
fsample = 1./dt;

% freq content
%figure(1)
%subplot(2,1,1)
avfft = mean(abs(fft(pltdata)),2);

figure; plot((0:round(nt/2))*fsample/nt,20*log10(avfft(1:round(nt/2)+1))  )
axis tight
xlabel('Frequency (Hz)')
ylabel('Power dB')
title('Amplitude spectrum')
