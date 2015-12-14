function [sample_p, sample_t, trace_p, trace_t] = peak_picker( trace )
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
% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
% github
% https://github.com/AnalysePrestackSeismic/
% ------------------ FUNCTION DEFINITION ---------------------------------
% Function Description
% this function picks peaks and troughs on seperate arrays on a trace
% 
%%

% find the first and second derivatives of the max
first_deriv = diff(trace);
second_deriv = diff(trace,2);

% apply a signum filter to get samples at zero crossings and make only 1
% and 0's
sign_1_deriv = sign(first_deriv);
sign_1_deriv(sign_1_deriv == 0) = 1;
sign_1_deriv(sign_1_deriv <0) = 0;

% find the point where sign of 1st deriv changes
diffsign = diff(sign_1_deriv);
mdiffsign = diffsign;

% set the point to zero where second derivative is positive and pad, this
% finds the peaks in the max dataset
diffsign(sign(second_deriv) > 0) = 0;
diffsign = [0;diffsign];

% set the point to zero where second derivative is positive and pad, this
% finds the mins in the max dataset
mdiffsign(sign(second_deriv) <= 0) = 0;
mdiffsign = [0;mdiffsign];

sampling_scale=1:length(trace);% Define the sampling axis
sampling_scale=sampling_scale';% Transpose the sampling axis
sample_p=sampling_scale(diffsign==-1);
sample_t=sampling_scale(mdiffsign==1);
trace_p=trace(diffsign==-1);
trace_t=trace(mdiffsign==1);
end

