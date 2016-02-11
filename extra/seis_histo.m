% program to load DIGI data and plot IG crossplots
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

function [] = seis_histo( segyfile)

int = struct;                                           % Intialize structure for intercept

int.filename=segyfile;       % Construct the file name , path for the intercept file

[int.meta, ~, int.traces]=segy_to_mat('189','193',int.filename);                                 % Load the intercept file for the block

figure;
hist(int.traces(:,1000),100)

% [hist_t,~,hist_axis_t,~]=histcn([[grad_trace_t_proj(ll_t);-grad_dp_lim;grad_dp_lim] [int_trace_t(ll_t);-int_dp_lim;int_dp_lim]],n_bin,n_bin); % Make a 2 D histogram for intercept  and gradient crossplot
% figure(5);
% pcolor(hist_axis_t{2},hist_axis_t{1},hist_t(1:n_bin,1:n_bin));                                                                                             % plot the 2D histogram against the bin centres
% hold on;
% plot(hist_axis_pts{2},(slope_z(ic)*hist_axis_pts{2}+grad0),'-k','LineWidth',3);                                                                              % Plot the trendline
% plot(hist_axis_pts{2},(slope_z_p(ic)*hist_axis_pts{2}+grad0_z_p(ic)),'-r','LineWidth',3);
% plot(hist_axis_pts{2},(slope_z_t(ic)*hist_axis_pts{2}+grad0_z_t(ic)),'-m','LineWidth',3);
% xlim([-150 150]);
% ylim([-300 300]);
% title('Cross plot from troughs');xlabel('Intercept');ylabel('Gradient');
% hold off;

end
