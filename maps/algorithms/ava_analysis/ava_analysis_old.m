function [] = ava_analysis_old( job_meta_path,block,int_dir,grad_dir,maxz,win_chi,inc_chi,int_dp_lim,grad_dp_lim,post_process_flag)
% UNTITLED Summary of this function goes here

% Broad Algorithm Architecture


%INPUTS:
%   job_meta_path : path of the job meta file
%   int_dir
%   grad_dir
%   maxz : Maximum z value in m for depth or ms for time
%   block: 
% % Use Default Computation parameters if .....
% win_chi = 100;       % Calculation moving window for chi calculation (in m or ms)
% inc_chi = 30; 
% int_dp_lim=300;     % Maximum absolute value of intercept for intercept - gradient crossplots
% grad_dp_lim=500;    % Maximum absolute value of gradient for intercept - gradient crossplots% Increament of moving window for chi calculation (in m or ms)
% post_process_flag=1; % Post process if 1 dont if 0

%OUTPUTS:
% Time/depth trends in form of a structure..(need to describe this)
%%

%   Detailed explanation goes here
close all;
%tic;
job_meta = load(job_meta_path);                         % Load job meta information


int = struct;
grad= struct;

int.filename=strcat(job_meta.output_dir,'digi_results/',int_dir,'/',int_dir,'_block_',block,'.segy');
grad.filename=strcat(job_meta.output_dir,'digi_results/',grad_dir,'/',grad_dir,'_block_',block,'.segy');


[int.meta int.ilxl_bytes int.traces]=segy_to_mat('189','193',int.filename);
[grad.meta grad.ilxl_bytes grad.traces]=segy_to_mat('189','193',grad.filename);

cpa_directory = strcat(job_meta.output_dir,'crossplot_analysis/',maxz,win_chi,inc_chi,int_dp_lim,grad_dp_lim,post_process_flag); % The name of the crossplot analysis directory

% Convert strings to number
block=str2double(block);
maxz=str2double(maxz);
win_chi=str2double(win_chi);
inc_chi=str2double(inc_chi);
int_dp_lim=str2double(int_dp_lim);
grad_dp_lim=str2double(grad_dp_lim);
post_process_flag=str2double(post_process_flag);

maxzout=maxz/(job_meta.s_rate/1000); % COnverting maximum z from m or ms into no of samples


%% Parameters
% Plot Parameters
crossplot=1;
if crossplot ==1
    n_bin=51;           % Number of bins in crossplot space (please declare as an odd number)
    n_bin_pt=21;        % Number of bins in crossplot space specific to crossplot of peaks an troughs only (please declare as an odd number)
end
plot_result=1;          % Toggle this to 1 if you want to disply final results for the block
%Emperical trend defination by y=mx+c type equation, used for plots only
chi_slope=-2;
chi_int=19;

%% -------Define the grid -------------------------------

ns_win_chi = win_chi/(job_meta.s_rate/1000);                            % Convert window length into number of samples
ns_inc_chi = floor(inc_chi/(job_meta.s_rate/1000));                     % Convert increment length into number of samples
st_win_chi = 1:ns_inc_chi:(floor(maxzout/ns_inc_chi)*ns_inc_chi);       % Array of Start of windows
end_win_chi = st_win_chi+ns_win_chi;                                    % Array of End of windows
end_win_chi(end_win_chi>maxzout)=maxzout;
nwin_chi=length(st_win_chi);                                            % TNumber of chi calculation windows

%% initialize directories and  variables

% Make a directory for storing results if it doesnot exist already
if exist(cpa_directory,'dir') == 0
    mkdir(cpa_directory);
end

if block == job_meta.liveblocks(1)
    job_meta.cpa_directory=cpa_directory;
    save(job_meta_path,'-struct','job_meta','-v7.3');
    
end


% Trend from all points
slope_z = zeros (1,nwin_chi);
chi_z = zeros (1,nwin_chi);
dist_z=zeros (1,nwin_chi);
int_mean_z = zeros (1,nwin_chi);
grad_mean_z = zeros (1,nwin_chi);
int_std_z = zeros (1,nwin_chi);
grad_std_z = zeros (1,nwin_chi);

% Trend from peaks
slope_z_p = zeros (1,nwin_chi);
chi_z_p = zeros (1,nwin_chi);
dist_z_p=zeros (1,nwin_chi);
int_mean_z_p = zeros (1,nwin_chi);
grad_mean_z_p = zeros (1,nwin_chi);
int_std_z_p = zeros (1,nwin_chi);
grad_std_z_p = zeros (1,nwin_chi);


% Trend from troughs
slope_z_t = zeros (1,nwin_chi);
chi_z_t = zeros (1,nwin_chi);
dist_z_t=zeros (1,nwin_chi);
int_mean_z_t = zeros (1,nwin_chi);
grad_mean_z_t = zeros (1,nwin_chi);
int_std_z_t = zeros (1,nwin_chi);
grad_std_z_t = zeros (1,nwin_chi);

intercept_temp = int.traces(:,(sum(int.traces,1)~=0|sum(grad.traces,1)~=0));      % Eleminate the blanc columns
gradient_temp =  grad.traces(:,(sum(grad.traces,1)~=0|sum(int.traces,1)~=0));     % Eleminate the blanc columns


angle=0;grad0=0;
angle_p=0;grad0_p=0;
angle_t=0;grad0_t=0;

z_axis=((st_win_chi+end_win_chi)/2)*(job_meta.s_rate/1000);                                         % Define the Z axis
chi_z_emperical=z_axis*chi_slope/1000 +chi_int ;   % Emperical chi trend used in DIGI start model and EER projection


%% Main computation Loop through the designed moving window
for ic=1:(nwin_chi-1)
  
    intercept_ic = intercept_temp(st_win_chi(ic):end_win_chi(ic),:);                            % Intercept values in the current window
    gradient_ic = gradient_temp(st_win_chi(ic):end_win_chi(ic),:);                              % Gradient values in the current window
    intercept_ic_uw_raw =intercept_ic(:);                                                       % Unwrap the intercept matrix
    gradient_ic_uw_raw = gradient_ic(:);                                                        % Unwrap the gradient matrix
    %fprintf('\n Window : %g',ic);
    if (max(intercept_ic_uw_raw)>0 && max(gradient_ic_uw_raw)>0)
        intercept_ic_uw =intercept_ic_uw_raw((intercept_ic_uw_raw~=0)|(gradient_ic_uw_raw~=0)); % Remove the points where both intercept and gradient are zero
        gradient_ic_uw = gradient_ic_uw_raw((intercept_ic_uw_raw~=0)|(gradient_ic_uw_raw~=0));  % Remove the points where both intercept and gradient are zero        
        clear intercept_ic_uw_raw gradient_ic_uw_raw intercept_ic gradient_ic;        
        ll=(abs(intercept_ic_uw)<int_dp_lim) & (abs(gradient_ic_uw)<grad_dp_lim);               %Logical for clipping intercept and gradient
        if(sum(ll)>1e4)
            [int_mean_z(ic) grad_mean_z(ic) int_std_z(ic) grad_std_z(ic) slope_z(ic) grad0 angle chi_z(ic) dist_z(ic)]=ava_statistics(intercept_ic_uw(ll),gradient_ic_uw(ll));     % Calculate ava stats from all points in intercept and gradient
        end
        
        if crossplot==1
            figure(1);            
            [hist_all,hist_axis]=hist3([[gradient_ic_uw(ll);-grad_dp_lim;grad_dp_lim] [intercept_ic_uw(ll);-int_dp_lim;int_dp_lim]],[n_bin n_bin]); % Make a 2 D histogram for intercept  and gradient crossplot
            hist_all=hist_all/max(max(hist_all));                                                                                                   % Normalize the histogram crossplot
            pcolor(hist_axis{2},hist_axis{1},hist_all);                                                                                             % plot the 2D histogram against the bin centres
            hold on;
            plot(hist_axis{2},(slope_z(ic)*hist_axis{2}+grad0),'--w');                                                                              % Plot the trendline
            hold off;
            xlabel('Intercept');ylabel('Gradient');                                                                                                 % Label the plot
            title(' Intercept-Gradient Crossplot (point density)');
        end
        
        %---------------Find Peaks and Troghs and do AVA staistics on them----------------
        
        [int_sample_p  int_sample_t int_trace_p int_trace_t]=peaks_and_troughs(intercept_ic_uw); % Find peaks and troughs in intercept trace
        [grad_sample_p  grad_sample_t grad_trace_p grad_trace_t]=peaks_and_troughs(gradient_ic_uw); % FInd peaks and troughs in gradient trace
        
        % Account for the fact that their are more peaks and troughs in the
        % intercept trace than the gradient trace and they might be at
        % slightly different locations as well
        
        %-----------------project gradient------------------        
        %         grad_trace_p_proj =  interp1(grad_sample_p,grad_trace_p,int_sample_p,'nearest'); % Project (nearest neighbour map) the gradient trace on the grid of the intercept trace (for peaks)
        %         grad_trace_t_proj =  interp1(grad_sample_t,grad_trace_t,int_sample_t,'nearest'); % Project (nearest neighbour map) the gradient trace on the grid of the intercept trace (for troughs)
        %         grad_trace_p=grad_trace_p_proj;
        %         grad_trace_t=grad_trace_t_proj;
        %                      
        %         % Extrapolate the ends to ensure edge effects and and ensure NaNs are removed
        %         grad_trace_p(end)=grad_trace_p(end-1);
        %         grad_trace_t(end)=grad_trace_t(end-1);
        %         grad_trace_p(1)=grad_trace_p(2);
        %         grad_trace_t(1)=grad_trace_t(2);

        %         grad_trace_p=grad_trace_p((~isnan(grad_trace_p)));   % remove nans from the gradient peak trace
        %         int_trace_p=int_trace_p((~isnan(grad_trace_p)));               % remove corresponding intercept points in the intercept peak trace
        %         grad_trace_t=grad_trace_t((~isnan(grad_trace_t)));   % remove nans from the gradient tough trace
        %         int_trace_t=int_trace_t((~isnan(grad_trace_t)));               % remove corresponding intercept points gradient tough trace
        %-------------------------------------------------------------


        %-----------------project intercept------------------
        int_trace_p_proj =  interp1(int_sample_p,int_trace_p,grad_sample_p,'nearest'); % Project (nearest neighbour map) the gradient trace on the grid of the intercept trace (for peaks)
        int_trace_t_proj =  interp1(int_sample_t,int_trace_t,grad_sample_t,'nearest'); % Project (nearest neighbour map) the gradient trace on the grid of the intercept trace (for troughs)
        int_trace_p=int_trace_p_proj;
        int_trace_t=int_trace_t_proj;
                     
        % Extrapolate the ends to ensure edge effects and and ensure NaNs are removed
        int_trace_p(end)=int_trace_p(end-1);
        int_trace_t(end)=int_trace_t(end-1);
        int_trace_p(1)=int_trace_p(2);
        int_trace_t(1)=int_trace_t(2);
        
        grad_trace_p=grad_trace_p((~isnan(int_trace_p)));   % remove nans from the gradient peak trace
        int_trace_p=int_trace_p((~isnan(int_trace_p)));               % remove corresponding intercept points in the intercept peak trace
        grad_trace_t=grad_trace_t((~isnan(int_trace_t)));   % remove nans from the gradient tough trace
        int_trace_t=int_trace_t((~isnan(int_trace_t)));               % remove corresponding intercept points gradient tough trace
        %-----------------------------------------------------
                       
        ll_p=(abs(int_trace_p)<int_dp_lim) & (abs(grad_trace_p)<grad_dp_lim); %Logical for clipping intercept and gradient
        ll_t=(abs(int_trace_t)<int_dp_lim) & (abs(grad_trace_t)<grad_dp_lim); %Logical for clipping intercept and gradient        
       % Do AVA statistics on the peaks
        if(sum(ll_p)>100)
            [int_mean_z_p(ic) grad_mean_z_p(ic) int_std_z_p(ic) grad_std_z_p(ic) slope_z_p(ic)  grad0_p angle_p chi_z_p(ic) dist_z_p(ic)]=ava_statistics(int_trace_p(ll_p),grad_trace_p(ll_p));
        end
        % Do ava statistics on troughs
        if(sum(ll_t)>100)
            [int_mean_z_t(ic) grad_mean_z_t(ic) int_std_z_t(ic) grad_std_z_t(ic) slope_z_t(ic) grad0_t angle_t chi_z_t(ic) dist_z_t(ic)]=ava_statistics(int_trace_t(ll_t),grad_trace_t(ll_t));
        end
        
                
        if crossplot==1
            figure(2);
                       
            hist_all_p=hist3([[int_trace_p(ll_p);-int_dp_lim;int_dp_lim]  [-1*grad_trace_p(ll_p);-grad_dp_lim;grad_dp_lim] ],[n_bin_pt n_bin_pt]);
            hist_all_p=hist_all_p';
            scatter(int_trace_p,grad_trace_p,2,[0 0 1]);xlim([-int_dp_lim int_dp_lim]);ylim([-grad_dp_lim grad_dp_lim]); % Plot peaks in blue
            %contour(hist_all_p,2);colormap jet;
            %imagesc(hist_all);
            
            hold on;
            
            hist_all_t=hist3([[int_trace_t(ll_t);-int_dp_lim;int_dp_lim]  [-1*grad_trace_t(ll_t);-grad_dp_lim;grad_dp_lim] ],[n_bin_pt n_bin_pt]);
            hist_all_t=hist_all_t';
            scatter(int_trace_t,grad_trace_t,2,[0 1 0]);xlim([-int_dp_lim int_dp_lim]);ylim([-grad_dp_lim grad_dp_lim]); % Plot peaks in green
            
            hist_all_pt=(hist_all_p+hist_all_t)/2;
            %contour(hist_all_t,2);colormap jet;
            %imagesc(hist_all_pt);
            
            hold on;
            % Plot trend lines through crossplots
            int_axis =[-int_dp_lim 0 int_dp_lim];
            plot(int_axis,(slope_z(ic)*int_axis+grad0),'-r');hold on;
            plot(int_axis,(slope_z_p(ic)*int_axis+grad0_p),'-b');hold on;
            plot(int_axis,(slope_z_t(ic)*int_axis+grad0_t),'-g');hold off;
            xlabel('Intercept');ylabel('Gradient');  title ('Intercept-Gradient scatter plot of peak and trough extremas');    
        end
    end
end
close all;
%toc;



%% Store the redults in a systematic fashion


clear intercept_ic gradient_ic intercept_ic_uw gradient_ic_uw intercept_temp gradient_temp ll ll_p ll_t crossplot;

% Find and store cetral inline and xline location

loc= job_meta.block_keys(block,:);
results_AVA.inl_central= floor((loc(1)+loc(2))/2);
results_AVA.xl_central = floor((loc(3)+loc(4))/2);

%results_AVA.header = 'intecept std deviation \n Gradient std deviation \n slope of trendline \n chi angle \n Emperical chi angle \n Distance of trend from origin';
results_AVA.header = [{'intecept_mean'},{'gradient_mean'},{'intecept_std_deviation'},{'gradient_std_deviation'},{'slope_of_trendline'},{'chi_angle_trend'},{'Dist_trend_from_origin'}];
results_AVA.zaxis=z_axis;
results_AVA.data = [int_mean_z;grad_mean_z;int_std_z;grad_std_z;slope_z;chi_z;dist_z];                     % Write various results into one matrix (for all points)
results_AVA.data_peaks=[int_mean_z_p;grad_mean_z_p;int_std_z_p;grad_std_z_p;slope_z_p;chi_z_p;dist_z_p];       % Write various results into one matrix (for peaks)
results_AVA.data_troughs=[int_mean_z_t;grad_mean_z_t;int_std_z_t;grad_std_z_t;slope_z_t;chi_z_t;dist_z_t];     % Write various results into one matrix (for troughs)

% Post process if required
if post_process_flag==1
    results_AVA.data=postprocess_results(results_AVA.data);
    results_AVA.data_peaks=postprocess_results(results_AVA.data_peaks);
    results_AVA.data_troughs=postprocess_results(results_AVA.data_troughs);
end

% Linear fit of chi trend to only peaks and troughs

chi_fit_p=polyfit(z_axis(~isnan(chi_z_p)),chi_z_p(~isnan(chi_z_p)),1);
chi_lin_p=chi_fit_p(1)*z_axis+chi_fit_p(2);

chi_fit_t=polyfit(z_axis(~isnan(chi_z_t)),chi_z_t(~isnan(chi_z_t)),1);
chi_lin_t=chi_fit_t(1)*z_axis+chi_fit_t(2);

% Store results
save(strcat(cpa_directory,'crossplot_analysis_block',num2str(block),'.mat'),'-struct','results_AVA','-v7.3');

%% --------------Plot Final Results --------------------------------------------------------------

if plot_result==1
    figure (3);
    
    for plot_inc=1:7
        subplot(1,7,plot_inc);
        plot(results_AVA.data(plot_inc,:),z_axis,'-r','LineWidth',3); set(gca,'YDir','rev'); % Plot depth trend of all points in red
        hold on;
        plot(results_AVA.data_peaks(plot_inc,:),z_axis,'-b','LineWidth',3); set(gca,'YDir','rev');% Plot depth trend of all peaks in blue
        plot(results_AVA.data_troughs(plot_inc,:),z_axis,'-g','LineWidth',3); set(gca,'YDir','rev');% Plot depth trend of all troughs in blue
        ylim([min(z_axis) max(z_axis)]);ylabel('depth/time');
        switch  plot_inc
            case 1
                title('Intercept Mean');
            case 2
                title('Gradient Mean');
                
            case 3
                title('Intercept Std dev');
            case 4
                title('Gradient Std dev');
            case 5
                title('Slope');
                xlim([-50 50]);
                
            case 6
                %plot(chi_lin_p,z_axis,'--b','LineWidth',3); set(gca,'YDir','rev');
                %plot(chi_lin_t,z_axis,'--g','LineWidth',3);
                plot(chi_z_emperical,z_axis,'--k','LineWidth',3);xlim([-40 40]);
                title('Chi Trend');
            case 7
                title('y-intercept of trend');xlim([-1000 1000]);
        end
        
    end
    set(gca,'YDir','rev');
    hold off;
end
end

%%####################Function to calculate ava statistics from intercept and gradient volume segments####################
function [int_mean_z grad_mean_z int_std_z grad_std_z slope_z grad0_z angle chi_trend dist_trend]=ava_statistics(intercept,gradient)
%% Funtion to do statistics on IG crossplot
% INPUT:
%   method = 1 for polyfit(y=mx +c) 2 for L2 for (y=mx) ype problems
%   intercept: Vector of intercept points
%   gradient: Vector of gradient points

% OUTPUT:

int_mean_z=mean(intercept);
grad_mean_z=mean(gradient);

%Parameters
plotQC=0;                                                           % Flag for making plots for QC on
fit_type = 2;                                                       % Define type of fitting to intercept Gradient Crossplot (values 1 , 2 or 3)
fraction_th=0.3;                                                    % Threshold to filterout points from normalized histogram
%---------------------
int_std_z  = std (intercept);                                       % Calculate intercept standard deviation
grad_std_z  = std (gradient);                                       % Calculate gradient standard deciation
%slope_z=lsqr(gradient,intercept);                                  % Calculating the least square fit
nbinfit=51;                                                         % Declare Number of bins for the 2D histogram of the intercept Gardient Crossplot
[hist_intgrad,~,hist_axis,~]=histcn([gradient intercept],nbinfit,nbinfit);         % Make a 2D histogram from the intercept gradient crossplot
%[hist_intgrad,hist_axis]=hist3([gradient intercept],[nbinfit nbinfit]);         % Make a 2D histogram from the intercept gradient crossplot
hist_intgrad=hist_intgrad/max(max(hist_intgrad));                   % Normalize the 2D histogramof IG crossplot

int_axis =hist_axis{2};
grad_axis=hist_axis{1};


%---------------Fit to IG crossplot--------------
% Fit to scatter directly
if(fit_type==1)
    intercept_w =intercept;
    gradient_w =gradient;
    clear intercept gradient;
end

% Fit to filtered scater, filter using histogram threshold
if(fit_type==2)
    int_inc=((max(intercept)-min(intercept))/nbinfit);                  % Calculate intercept axis increment
    grad_inc=((max(gradient)-min(gradient))/nbinfit);                   % Calculate gradient axis increment
    int_axis = min(intercept)+(0:(nbinfit-1)+0.5)*int_inc;              % Design  the intercept axis
    grad_axis = min(gradient)+(0:(nbinfit-1)+0.5)*grad_inc;             % Design the gradient axis
    ind_int=floor((intercept-min(intercept))/int_inc);
    ind_grad=floor((gradient-min(gradient))/grad_inc);
    ind_int(ind_int<1)=1;ind_int(ind_int>nbinfit)=nbinfit;
    ind_grad(ind_grad<1)=1;ind_grad(ind_grad>nbinfit)=nbinfit;
    
    weight_ig=ones(1,length(intercept));
    for f1=1:length(ind_grad)
        weight_ig(f1)=hist_intgrad(ind_grad(f1),ind_int(f1));
    end
    weight_ig=weight_ig';
    intercept_w=intercept(weight_ig>fraction_th);
    gradient_w=gradient(weight_ig>fraction_th);
end

% Fit to histogram bin centres weigthed by the histogram value in the fit
if(fit_type==3)
    k=0.05;
    n_ig_w=length(int_axis)*length(grad_axis)*ceil(1/k); % Maximum possible number of reconstructed I and G
    intercept_w=ones(n_ig_w,1)*999.99;
    gradient_w=ones(n_ig_w,1)*999.99;
    c1=1;
    for f1=1:length(int_axis)
        for f2=1:length(grad_axis)
            if hist_intgrad(f2,f1)>fraction_th
                for fr=k:k:hist_intgrad(f2,f1)
                    intercept_w(c1)=int_axis(f1);
                    gradient_w(c1)=grad_axis(f2);
                    c1=c1+1;
                end
            end
        end
    end
    intercept_w=intercept_w(intercept_w~=999.99);
    gradient_w =gradient_w(gradient_w~=999.99);
end
%--------------------------------------------------

p1=polyfit(intercept_w,gradient_w,1);% Solve minimizing  y-error
%p=polyfit(gradient,intercept,1);
slope1=p1(1);
g_int1=p1(2);

p2=polyfit(gradient_w,intercept_w,1);% Solve minimizing  x-error
slope2=1/p2(1);
g_int2=-p2(2)/p2(1);

% Average the two fits
slope_z=(slope1+slope2)/2;
grad0_z=(g_int1+g_int2)/2;

%if(slope_z
%--------------Plots for QC -------------------------------------
if plotQC==1
    figure(11);pcolor(int_axis,grad_axis,hist_intgrad);
    hold on;plot(intercept_w,(slope_z*intercept_w+grad0_z),'--w');hold off;
    
    figure(12);
    scatter(intercept_w,gradient_w,1);
    xlim([min(intercept) max(intercept)]);
    ylim([min(gradient) max(gradient)]);
    hold on;
    plot(intercept_w,(slope_z*intercept_w+grad0_z),'--r');
    hold off;
end
%--------------------------------------------------------------   

angle= atand(slope_z);
% Calculating the angle of slope least square fit
%chi_trend =90+angle;
if angle>0
    chi_trend=angle-90;
else
    chi_trend=90+angle;
end
% if chi_trend >90
%    chi_trend = -(chi_trend-90);
% end
%chi_trend =(90+angle)*(angle<0)+(-90+angle)*(angle>0);                                                         % Transform into chi angle
%chi_trend = -1*angle;                                                         % Transform into chi angle
%dist_trend=p(2)*cosd(chi_trend);
dist_trend=grad0_z;%*cosd(angle);
end

%% ########################### Function to return the peaks and troughs in a trace and their location as
% samples
function [sample_p sample_t trace_p trace_t]=peaks_and_troughs(trace)
% find the max of the data across the gather and smooth
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

% code to make a trace into a blocky trace honouring peaks and troughs only
% trace_p_blocky=interp1(sample_p,trace_p,sampling_scale,'nearest');
% trace_t_blocky=interp1(sample_t,trace_t,sampling_scale,'nearest');
% trace_blocky =interp1([sample_p sample_t],[trace_p trace_t],sampling_scale,'nearest');
end
%% ##################function to postprocess poorly conditioned data and replace it by median value of the result array
function[M2]= postprocess_results(M)
flag_baddata=(M(3,:)<0.1*max(M(3,:)))|(M(4,:)<0.1*max(M(4,:)));
M2=M;
for p=5:7
    med_val=median(M(p,:));
    M2(p,flag_baddata)=med_val;    
end
end
