function [vint] = vrms_vint_dix_js(z_in,vrms_in,weight)
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
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

z_step = 4; % 4 ms is 1 index step

% Model is Vint^2 at denser sampling

% Grid onto discrete index locations
%min_z_in = min(z);
%max_z_in = max(z);
z_in = 1+(z_in/z_step);
min_z_out = 0;
max_z_out = 1240;

z_out = min_z_out:1:max_z_out;
n_out = length(z_out);

z_in_round = round(z_in);
z_in_round(end) = max_z_out;
n_in = length(z_in);
vrms_in_interp = interp1(z_in,vrms_in,z_in_round,'spline');

data = z_in_round.*(vrms_in_interp.^2); % Data is t*Vrms^2 where k is sample point

for iz = 1:n_in
    C(iz,1:z_in_round(iz)) = ones(1,z_in_round(iz));
end

%D = spdiags([weight*ones(max_z_out,1) -2*weight*ones(max_z_out,1) weight*ones(max_z_out,1)],[-1 0 1],max_z_out,max_z_out);
D = spdiags([weight*-1*ones(max_z_out,1) weight*ones(max_z_out,1)],[0 1],max_z_out,max_z_out);

vint = real(sqrt(lsqr([C;D],[data;zeros(max_z_out,1)])));
%vint = real(sqrt(lsqr([C],[data])));
%plot(sqrt(m),-1*(0:1:max_z_out-1))

end

