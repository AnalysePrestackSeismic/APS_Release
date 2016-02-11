function [ freq_pairs ] = spec_cos_taper( f1,f2,f3,f4 )
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
%generate freq, amplitude pairs from 4 corner points
%using supplied 4 corner frequencies
%note that this is hardwired for 1 Hz intervals which may not be good
%enough for low frequencies

freq_pairs = zeros(f4,2); % initialise output

freq_pairs(:,1) = [1:f4]; % set frequency values to range of corner points

% frequencies between f1 and f2 set to cosine taper up from 0 to 1
% frequencies between f2 and f3 are set to 1
% frequencies between f3 and f4 set to cosine taper down from 1 to 0

freq_pairs(f1:f2,2) = fliplr(0.5*(1+cos(-pi/(f2-f1).*([f1:f2]-f1))));
freq_pairs(f2:f3,2) = 1;
freq_pairs(f3:f4,2) = 0.5*(1+cos(-pi/(f4-f3).*([f3:f4]-f3)));

% plot the output

figure; plot(freq_pairs(:,1),freq_pairs(:,2))


end

