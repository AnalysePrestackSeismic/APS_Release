function [] = ava_analysis_aperture( job_meta_path_intercept,job_meta_path_gradient,block,maxz,win_chi,inc_chi,varargin)
% UNTITLED Summary of this function goes here

% Broad Algorithm Architecture


%INPUTS:
%   job_meta_path : path of the job meta file
%   int_dir
%   grad_dir
%   maxz : Maximum z value in m for depth or ms for time
%   block:
%   Use Default Computation parameters if .....
%   win_chi = 100;       % Calculation moving window for chi calculation (in m or ms)
%   inc_chi = 30;

%   varagin{1} =grad_dp_lim varagin{2} =int_dp_lim, varagin{3} =grad_dp_lim 

% aperture

%   int_dp_lim=300;     % Maximum absolute value of intercept for intercept -
%   gradient crossplots  .  Initialize this a 0 if you want the computer to
%   decide this on every window (recommemded)

%   grad_dp_lim=500;    % Maximum absolute value of gradient for intercept -
%   gradient crossplots% Increament of moving window for chi calculation (in
%   m or ms). .  Initialize this a 0 if you want the computer to
%   decide this on every window (recommemded) 



%OUTPUTS:
%   Time/depth trends in form of a structure..(need to describe this)
%% Handle User Defined Parameters

% Detailed explanation goes here
close all;
%tic;
job_meta_path = job_meta_path_intercept;                % Hoping that the intercept and gradient will have same geometry 
job_meta = load(job_meta_path);                         % Load job meta information

job_meta_int = load(job_meta_path_intercept);
job_meta_grad =load(job_meta_path_gradient);
int = struct;                                           % Intialize structure for intercept
grad= struct;                                           % Intialize structure for gradient

[~, int.traces, int.ilxl_read] = node_segy_read(job_meta_path_intercept,'1',block);
[~, grad.traces, grad.ilxl_read] = node_segy_read(job_meta_path_gradient,'1',block);

% int.filename=strcat(job_meta.output_dir,'digi_results/',int_dir,'/',int_dir,'_block_',block,'.segy');       % Construct the file name , path for the intercept file
% grad.filename=strcat(job_meta.output_dir,'digi_results/',grad_dir,'/',grad_dir,'_block_',block,'.segy');    % Construct the file name , path for the igradient file
% 

% If the advanced parameters are defined by user extract it from varargin
% else initialize them as zero
if ~isempty(varargin)
    aperture =str2double(varargin{1});
    
    if length(varagin)>1
        int_dp_lim =str2double(varargin{2});
        grad_dp_lim =str2double(varargin{3});
    else
        int_dp_lim = 0;
        grad_dp_lim = 0;
    end
    
else
    aperture =0;                                                                                            %Default aperture is zero ( => Neighbourhood Radius of zero)
    int_dp_lim = 0;
    grad_dp_lim = 0;
end

% [int.meta int.ilxl_bytes int.traces]=segy_to_mat('189','193',int.filename);                                 % Load the intercept file for the block
% [grad.meta grad.ilxl_bytes grad.traces]=segy_to_mat('189','193',grad.filename);                             % Load the gradient file for the block

cpa_directory = strcat(job_meta.output_dir,'crossplot_analysis','_','_',num2str(maxz),'_',num2str(win_chi),'_',num2str(inc_chi),'/'); % The name of the crossplot analysis directory

% Convert strings to number
block=str2double(block);
maxz=str2double(maxz);
win_chi=str2double(win_chi);
inc_chi=str2double(inc_chi);
int_dp_lim=str2double(int_dp_lim);
grad_dp_lim=str2double(grad_dp_lim);

maxzout=maxz/(job_meta.s_rate/1000); % Converting maximum z from m or ms into no of samples
if  maxzout==0 || maxzout> size(int.traces,1)
    
    maxzout = size(int.traces,1);
end

%% Hardwired Parameters
% Plot Parameters
crossplot=1;
qcplot=0;
if crossplot ==1
    n_bin=51;           % Number of bins in crossplot space (please declare as an odd number)
end
plot_result=1;          % Toggle this to 1 if you want to disply final results for the block


%Emperical trend defination by y=mx+c type equation, used for plots only
chi_slope=-2;
chi_int=19;

percentile_keep = 90; % This determines the percentile of intercept and gradient windows to keep. Calculated for every window
padding = 50;       % water bottom pick padding

%% -------Define the grid -------------------------------
% % Pick water bottom or use a pre picked water bottom horizon
% if isfield(job_meta, 'wb_path')
%     wb_idx_in = dlmread(job_meta.wb_path);
%     % col 1 inline
%     % col 2 xline
%     % col 3 twt
%     if job_meta.is_gather == 1
%         [~,locations] = ismember(results_out{1,2}{1,1}(1:length(offset):end,:),wb_idx_in(:,1:2),'rows');
%     else
%         error('should be gathers');
%         %[~,locations] = ismember(ilxl_read{1}(1:end,:),wb_idx_in(:,1:2),'rows');
%     end
%     %wb_idx = zeros(size(traces{vol_index_wb},2),1);
%     zero_loc = locations ~= 0;
%     % make a zeros for each trace
%     xi = (1:n_traces)';
%     x = xi(zero_loc)';
%     wb_idx = interp1(x,wb_idx_in(locations(zero_loc),3),xi);
%     clear wb_idx_in
%     
%     wb_idx = (wb_idx./(job_meta.s_rate/1000))';
%     
%     wb_idx = round(wb_idx-padding);
%     wb_idx(isnan(wb_idx)) = 1;
%     wb_idx(wb_idx < 1) = 1;
%     %win_sub = bsxfun(@plus,wb_idx,(0:job_meta.n_samples{1}-max(wb_idx))');
%     
%     %win_ind = bsxfun(@plus,win_sub,(0:job_meta.n_samples{1}:job_meta.n_samples{1}*(size(traces{vol_index_wb},2)-1)));
% else
%    [wb_idx] = water_bottom_picker(stackin_freq,padding);
    [wb_idx] = water_bottom_picker(int.traces,padding);
    wb_idx(wb_idx < 0) = 1;
    wb_idx = wb_idx+65;
% end
%--------------------------------------------------------------------
% wb_sample= ceil(job_meta.wb_z_avg(block));
wb_sample =wb_idx ;
ns_win_chi = win_chi/(job_meta.s_rate/1000);                            % Convert window length into number of samples
ns_inc_chi = floor(inc_chi/(job_meta.s_rate/1000));                     % Convert increment length into number of samples
st_win_chi = 1:ns_inc_chi:(floor(maxzout/ns_inc_chi)*ns_inc_chi);       % Array of Start of windows
end_win_chi = st_win_chi+ns_win_chi;                                    % Array of End of windows
end_win_chi(end_win_chi>maxzout)=maxzout;
nwin_chi=length(st_win_chi);                                            % Number of chi-calculation windows

%st_win_chi(st_win_chi<wb_sample)=wb_sample;                             % Limit start of windows at the water bottom
%end_win_chi(end_win_chi<wb_sample)=wb_sample;                           % Limit end of windows at the water bottom

%% initialize directories and  variables and condition other inputs

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

% Trend from extremas
slope_z_pts = zeros (1,nwin_chi);
chi_z_pts = zeros (1,nwin_chi);
dist_z_pts=zeros (1,nwin_chi);
int_mean_z_pts = zeros (1,nwin_chi);
grad_mean_z_pts = zeros (1,nwin_chi);
int_std_z_pts = zeros (1,nwin_chi);
grad_std_z_pts = zeros (1,nwin_chi);

% Trend from peaks
slope_z_p = zeros (1,nwin_chi);
chi_z_p = zeros (1,nwin_chi);
dist_z_p=zeros (1,nwin_chi);
int_mean_z_p = zeros (1,nwin_chi);
grad_mean_z_p = zeros (1,nwin_chi);
int_std_z_p = zeros (1,nwin_chi);
grad_std_z_p = zeros (1,nwin_chi);
grad0_z_p=zeros (1,nwin_chi);


% Trend from troughs
slope_z_t = zeros (1,nwin_chi);
chi_z_t = zeros (1,nwin_chi);
dist_z_t=zeros (1,nwin_chi);
int_mean_z_t = zeros (1,nwin_chi);
grad_mean_z_t = zeros (1,nwin_chi);
int_std_z_t = zeros (1,nwin_chi);
grad_std_z_t = zeros (1,nwin_chi);
grad0_z_t=zeros (1,nwin_chi);


intercept_temp = int.traces(:,(sum(int.traces,1)~=0|sum(grad.traces,1)~=0));      % Eliminate the blanc columns
intercept_temp(1:wb_sample+50,:)=0;                                               % Mute everything till the waterbottom
gradient_temp =  grad.traces(:,(sum(grad.traces,1)~=0|sum(int.traces,1)~=0));     % Eliminate the blanc columns
gradient_temp(1:wb_sample+50,:)=0;                                                % Mute everything till the waterbottom

grad0=0;
grad0_p=0;
z_axis=((st_win_chi+end_win_chi)/2)*(job_meta.s_rate/1000);                       % Define the Z axis
chi_z_emperical=z_axis*chi_slope/1000 +chi_int ;   % Emperical chi trend used in DIGI start model and EER projection
flag_lim=0;
int_dp_lim_bk=int_dp_lim;                           % Store defult value of intercept limit
grad_dp_lim_bk=grad_dp_lim;                         % Store default value of gradient limit

clear int grad;

flag_check=ones((nwin_chi),1);
% nsamples_th=win_chi*1e2;                                                                            % Threshold number of samples to do computation for that window
nsamples_th=win_chi*2;   
%% Main computation Loop through the designed moving window


for ic=1:(nwin_chi-1)
    intercept_ic = intercept_temp(st_win_chi(ic):end_win_chi(ic),:);                                % Intercept values in the current window
    gradient_ic = gradient_temp(st_win_chi(ic):end_win_chi(ic),:);                                  % Gradient values in the current window
    intercept_ic_uw_raw =intercept_ic(:);                                                           % Unwrap the intercept matrix
    gradient_ic_uw_raw = gradient_ic(:);                                                            % Unwrap the gradient matrix
    
    if (max(abs(intercept_ic_uw_raw))>0 && max(abs(gradient_ic_uw_raw))>0)
        % Automate intercept and gradient limits
        
        intercept_ic_uw =intercept_ic_uw_raw((intercept_ic_uw_raw~=0)|(gradient_ic_uw_raw~=0));     % Remove the points where both intercept and gradient are zero
        gradient_ic_uw = gradient_ic_uw_raw((intercept_ic_uw_raw~=0)|(gradient_ic_uw_raw~=0));      % Remove the points where both intercept and gradient are zero
        
        % ###################FIND A SENSIBLE RANGE FOR IG CROSSPLOT###########################
        if length(intercept_ic_uw)>nsamples_th
            %if (int_dp_lim==0||grad_dp_lim==0)
                flag_lim=flag_lim+1;                                                                                % Flag that the limits are being automated
                % Method : 1
                %             scale_lim=4;                                                                                % Scale factor for standard deviation into range
                %             [~,~,int_std_temp,grad_std_temp,~,~,~,~,~]=ava_statistics(intercept_ic_uw,gradient_ic_uw);  % Find Ava statistics on raw data and find standar deviation of intercept and gradient
                %             int_dp_lim = scale_lim*int_std_temp;                                                        % Scale standard deviation to make a range for intercept
                %             grad_dp_lim = scale_lim*grad_std_temp;                                                      % Scale standard deviation to make a range for gradient
                %             clear int_std_temp  grad_std_temp;
                
                % Method : 2
                %               scale_lim=0.1;
                %               int_dp_lim =scale_lim* max(abs(intercept_ic_uw));
                %               grad_dp_lim =scale_lim*max(abs(gradient_ic_uw));
                
                
                % Method: 3
                
                intercept_ic_uw_sort = sort(abs(intercept_ic_uw));                                                          % Sort the absolute values of intercept
                int_dp_lim2 = intercept_ic_uw_sort(floor((percentile_keep/100)*length(intercept_ic_uw_sort)));               % Find what should be the range to keep the right percentile of the points
                if qcplot==1
                    figure(21);
                    subplot(1,2,1)
                    plot(intercept_ic_uw_sort);title('Intercept');
                end
                clear intercept_ic_uw_sort;
                
                gradient_ic_uw_sort = sort(abs(gradient_ic_uw));
                grad_dp_lim2 = gradient_ic_uw_sort(floor((percentile_keep/100)*length(gradient_ic_uw_sort)));
                if qcplot==1
                    subplot(1,2,2)
                    plot(gradient_ic_uw_sort,'-r');title('Gradient');
                end
                clear gradient_ic_uw_sort;
                
                if((int_dp_lim2<int_dp_lim_bk)||int_dp_lim_bk==0)
                    int_dp_lim = int_dp_lim2;
                else
                    int_dp_lim = int_dp_lim_bk;
                end
                
                if((grad_dp_lim2<grad_dp_lim_bk)||grad_dp_lim_bk==0)
                    grad_dp_lim = grad_dp_lim2;
                else
                    grad_dp_lim = grad_dp_lim_bk;
                end
        end
        
        %############################### DO AVA STATISTICS #####################################
        clear intercept_ic_uw_raw gradient_ic_uw_raw intercept_ic gradient_ic;
        ll=(abs(intercept_ic_uw)<int_dp_lim) & (abs(gradient_ic_uw)<grad_dp_lim);                   %Logical for clipping intercept and gradient
        % Execute this point if you have at least threshold number of points points for the crossplot
        if(sum(ll)>(percentile_keep/100)*nsamples_th)
            [int_mean_z(ic),grad_mean_z(ic),int_std_z(ic),grad_std_z(ic),slope_z(ic),~,~, chi_z(ic), dist_z(ic)]=ava_statistics(intercept_ic_uw(ll),gradient_ic_uw(ll),1);        % Calculate ava stats from all points in intercept and gradient
            if crossplot==1
                figure(1);
                [hist_all,~,hist_axis,~]=histcn([[gradient_ic_uw(ll);-grad_dp_lim;grad_dp_lim] [intercept_ic_uw(ll);-int_dp_lim;int_dp_lim]],n_bin,n_bin);                      % Make a 2 D histogram for intercept  and gradient crossplot
                hist_all=hist_all/max(max(hist_all));                                                                                                                           % Normalize the histogram crossplot
                pcolor(hist_axis{2},hist_axis{1},hist_all(1:n_bin,1:n_bin));                                                                                                    % plot the 2D histogram against the bin centres
                hold on;
                plot(hist_axis{2},(slope_z(ic)*hist_axis{2}+grad0),'--w');                                                                                                      % Plot the trendline
                hold off;
                xlabel('Intercept');ylabel('Gradient');                                                                                                                         % Label the plot
                title(' Intercept-Gradient Crossplot (point density)');
            end
            
            %---------------Find Peaks and Troghs and do AVA staistics on them----------------
            
            [int_sample_p int_sample_t int_trace_p int_trace_t]=peaks_and_troughs(intercept_ic_uw);         % Find peaks and troughs in intercept trace
            %[grad_sample_p grad_sample_t grad_trace_p grad_trace_t]=peaks_and_troughs(gradient_ic_uw);     % Find peaks and troughs in gradient trace
            
            %-----------------project intercept------------------
            grad_trace_p_proj =gradient_ic_uw(int_sample_p);                                                % Project (nearest neighbour map) the gradient trace on the grid of the intercept trace (for peaks)
            grad_trace_t_proj =gradient_ic_uw(int_sample_t);                                                % Project (nearest neighbour map) the gradient trace on the grid of the intercept trace (for troughs)
            int_trace_pts=[int_trace_p;int_trace_t];
            grad_trace_pts=[grad_trace_p_proj;grad_trace_t_proj];
            
            nans_pts=isnan(grad_trace_pts)|isnan(int_trace_pts);
            grad_trace_pts=grad_trace_pts(~nans_pts);                                                       % remove nans from the gradient peak trace
            int_trace_pts=int_trace_pts(~nans_pts);                                                         % remove corresponding intercept points in the intercept peak trace
            ll_pts=(abs(int_trace_pts)<int_dp_lim) & (abs(grad_trace_pts)<grad_dp_lim);                     %Logical for clipping intercept and gradient
            
            % ------------Do AVA statistics on the peaks and troughs------------------------
            if(sum(ll_pts)>(percentile_keep/100)*nsamples_th*0.1)
                [int_mean_z_pts(ic),grad_mean_z_pts(ic),int_std_z_pts(ic),grad_std_z_pts(ic),slope_z_pts(ic),~,~, chi_z_pts(ic),dist_z_pts(ic)]=ava_statistics(int_trace_pts(ll_pts),grad_trace_pts(ll_pts),1);
                flag_check(ic)=0;
                
                % Plot progresss ofbackground chi computation
                if qcplot==1
                    figure(6);plot(chi_z);
                    hold on; plot(chi_z_pts,'-r');
                    hold off;
                end
                
            end
            
            % -------------Do AVA Statistics for Peaks-----------------------
            
            ll_p=(abs(int_trace_p)<int_dp_lim) & (abs(grad_trace_p_proj)<grad_dp_lim); %Logical for clipping intercept and gradient
            
            if(sum(ll_p)>(percentile_keep/100)*nsamples_th*0.03)
                [int_mean_z_p(ic),grad_mean_z_p(ic),int_std_z_p(ic),grad_std_z_p(ic),slope_z_p(ic),grad0_z_p(ic),~, chi_z_p(ic),dist_z_p(ic)]=ava_statistics(int_trace_p(ll_p),grad_trace_p_proj(ll_p),2);
                flag_check(ic)=0;
                
                % Plot progresss ofbackground chi computation
                if qcplot==1
                    figure(6);plot(chi_z);
                    hold on; plot(chi_z_p,'-g');
                    hold off;
                end
                
            end
            
            % -------------Do AVA Statistics for Peaks-----------------------
            
            ll_t=(abs(int_trace_t)<int_dp_lim) & (abs(grad_trace_t_proj)<grad_dp_lim); %Logical for clipping intercept and gradient
            
            if(sum(ll_p)>(percentile_keep/100)*nsamples_th*0.03)
                [int_mean_z_t(ic),grad_mean_z_t(ic),int_std_z_t(ic),grad_std_z_t(ic),slope_z_t(ic),grad0_z_t(ic),~, chi_z_t(ic),dist_z_t(ic)]=ava_statistics(int_trace_t(ll_t),grad_trace_t_proj(ll_t),2);
                flag_check(ic)=0;
                
                % Plot progresss ofbackground chi computation
                if qcplot==1
                    figure(6);plot(chi_z);
                    hold on; plot(chi_z_t,'-m');
                    hold off;
                end
                
            end
            
            
            %------------Display Crossplot-------------------------
            if crossplot==1
%                 figure(2);
%                 scatter(int_trace_pts,grad_trace_pts,2,[0 0 1]);xlim([-int_dp_lim int_dp_lim]);ylim([-grad_dp_lim grad_dp_lim]); % Plot peaks and troughs in blue
%                 hold on;
%                 % Plot trend lines through crossplots
%                 int_axis =[-int_dp_lim 0 int_dp_lim];
%                 plot(int_axis,(slope_z(ic)*int_axis+grad0),'-r');hold on;
%                 plot(int_axis,(slope_z_pts(ic)*int_axis+grad0_p),'-b');
%                 plot(int_axis,(slope_z_(ic)*int_axis+grad0_p),'-g');hold off;                
%                 xlabel('Intercept');ylabel('Gradient');  title ('Intercept-Gradient scatter plot of peak and trough extremas');
                
                
                figure(3);                
                [hist_pts,~,hist_axis_pts,~]=histcn([[grad_trace_pts(ll_pts);-grad_dp_lim;grad_dp_lim] [int_trace_pts(ll_pts);-int_dp_lim;int_dp_lim]],n_bin,n_bin); % Make a 2 D histogram for intercept  and gradient crossplot
                hist_pts=hist_pts/max(max(hist_pts));                                                                                                   % Normalize the histogram crossplot
                pcolor(hist_axis_pts{2},hist_axis_pts{1},hist_pts(1:n_bin,1:n_bin));                                                                                             % plot the 2D histogram against the bin centres
                hold on;
                plot(hist_axis_pts{2},(slope_z(ic)*hist_axis_pts{2}+grad0),'-w','LineWidth',3);                                                                              % Plot the trendline
                plot(hist_axis_pts{2},(slope_z_p(ic)*hist_axis_pts{2}+grad0_z_p(ic)),'-r','LineWidth',3); 
                plot(hist_axis_pts{2},(slope_z_t(ic)*hist_axis_pts{2}+grad0_z_t(ic)),'-k','LineWidth',3); 
                hold off;
                
                xlabel('Intercept');ylabel('Gradient');                                                                                                 % Label the plot
                title(' Intercept pts-Gradient pts Crossplot (point density)');hold off;
                
                
                
                [hist_p,~,hist_axis_p,~]=histcn([[grad_trace_p_proj(ll_p);-grad_dp_lim;grad_dp_lim] [int_trace_p(ll_p);-int_dp_lim;int_dp_lim]],n_bin,n_bin); % Make a 2 D histogram for intercept  and gradient crossplot
                figure(4);
                pcolor(hist_axis_p{2},hist_axis_p{1},hist_p(1:n_bin,1:n_bin));                                                                                             % plot the 2D histogram against the bin centres
                hold on;
                plot(hist_axis_pts{2},(slope_z(ic)*hist_axis_pts{2}+grad0),'-w','LineWidth',3);                                                                              % Plot the trendline
                plot(hist_axis_pts{2},(slope_z_p(ic)*hist_axis_pts{2}+grad0_z_p(ic)),'-r','LineWidth',3); 
                plot(hist_axis_pts{2},(slope_z_t(ic)*hist_axis_pts{2}+grad0_z_t(ic)),'-k','LineWidth',3); 
                title('Cross plot from peaks');xlabel('Intercept');ylabel('Gradient');   
                hold off;
                
                
                [hist_t,~,hist_axis_t,~]=histcn([[grad_trace_t_proj(ll_t);-grad_dp_lim;grad_dp_lim] [int_trace_t(ll_t);-int_dp_lim;int_dp_lim]],n_bin,n_bin); % Make a 2 D histogram for intercept  and gradient crossplot
                figure(5);
                pcolor(hist_axis_t{2},hist_axis_t{1},hist_t(1:n_bin,1:n_bin));                                                                                             % plot the 2D histogram against the bin centres
                hold on;
                plot(hist_axis_pts{2},(slope_z(ic)*hist_axis_pts{2}+grad0),'-w','LineWidth',3);                                                                              % Plot the trendline
                plot(hist_axis_pts{2},(slope_z_p(ic)*hist_axis_pts{2}+grad0_z_p(ic)),'-r','LineWidth',3); 
                plot(hist_axis_pts{2},(slope_z_t(ic)*hist_axis_pts{2}+grad0_z_t(ic)),'-k','LineWidth',3); 
                title('Cross plot from troughs');xlabel('Intercept');ylabel('Gradient');   
                hold off;
                
            end
        end
        if(flag_lim==1)
            int_dp_lim=0;
            grad_dp_lim=0;
        end
    end
    %end
    
end

close all;


%% Store the redults in a systematic fashion


clear intercept_ic gradient_ic intercept_ic_uw gradient_ic_uw intercept_temp gradient_temp ll ll_p ll_t crossplot int grad;

% Find and store cetral inline and xline location

loc= job_meta.block_keys(block,:);
results_AVA.inl_central= floor((loc(1)+loc(2))/2);
results_AVA.xl_central = floor((loc(3)+loc(4))/2);

%results_AVA.header = 'intecept std deviation \n Gradient std deviation \n slope of trendline \n chi angle \n Emperical chi angle \n Distance of trend from origin';
results_AVA.header = [{'intecept_mean'},{'gradient_mean'},{'intecept_std_deviation'},{'gradient_std_deviation'},{'slope_of_trendline'},{'chi_angle_trend'},{'Dist_trend_from_origin'}];
results_AVA.zaxis=z_axis;
results_AVA.data = [int_mean_z;grad_mean_z;int_std_z;grad_std_z;slope_z;chi_z;dist_z];                                      % Write various results into one matrix (for all points)
results_AVA.data_pts=[int_mean_z_pts;grad_mean_z_pts;int_std_z_pts;grad_std_z_pts;slope_z_pts;chi_z_pts;dist_z_pts];        % Write various results into one matrix (for peaks)
results_AVA.data_p=[int_mean_z_p;grad_mean_z_p;int_std_z_p;grad_std_z_p;slope_z_p;chi_z_p;dist_z_p];        % Write various results into one matrix (for peaks)
results_AVA.data_t=[int_mean_z_t;grad_mean_z_t;int_std_z_t;grad_std_z_t;slope_z_t;chi_z_t;dist_z_t];        % Write various results into one matrix (for peaks)



% Post process if required
% if post_process_flag==1
%     results_AVA.data=postprocess_results(results_AVA.data);
%     results_AVA.data_pts=postprocess_results(results_AVA.data_pts);
% end

% Linear fit of chi trend to only peaks and troughs

chi_fit_pts=polyfit(z_axis(~isnan(chi_z_pts)),chi_z_pts(~isnan(chi_z_pts)),1);
chi_lin_pts=chi_fit_pts(1)*z_axis+chi_fit_pts(2);

% Store results
save(strcat(cpa_directory,'crossplot_analysis_block',num2str(block),'.mat'),'-struct','results_AVA','-v7.3');

%% --------------Plot Final Results --------------------------------------------------------------

if plot_result==1
    figure (11);
    
    for plot_inc=1:7
        subplot(1,7,plot_inc);
        plot(results_AVA.data(plot_inc,:),z_axis,'-r','LineWidth',3); set(gca,'YDir','rev');                                % Plot depth trend of all points in red
        hold on;
        plot(results_AVA.data_pts(plot_inc,:),z_axis,'-b','LineWidth',3); set(gca,'YDir','rev');                            % Plot depth trend of all peaks and troughs in blue
        hold on;
        plot(results_AVA.data_p(plot_inc,:),z_axis,'-g','LineWidth',3); set(gca,'YDir','rev');                            % Plot depth trend of all peaks in green
        hold on;
        plot(results_AVA.data_t(plot_inc,:),z_axis,'-m','LineWidth',3); set(gca,'YDir','rev');                            % Plot depth trend of all peaks in magenta
        
        
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
                
                plot(chi_z_emperical,z_axis,'--k','LineWidth',3);
                title('Chi Trend');
            case 7
                %plot(flag_check,z_axis,'--k','LineWidth',3);
                title('distance of trend');%xlim([-2 2]);
        end
        
    end
    set(gca,'YDir','rev');
    hold off;
end
end

%% ####################Function to calculate ava statistics from intercept and gradient volume segments####################
function [int_mean_z,grad_mean_z,int_std_z,grad_std_z,slope_z,grad0_z,angle,chi_trend,dist_trend]=ava_statistics(intercept,gradient,method)
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
    ind_int=floor((intercept-min(intercept))/int_inc);                  % Create index/axis for intercept
    ind_grad=floor((gradient-min(gradient))/grad_inc);                  % Create index/axis for gradient
    ind_int(ind_int<1)=1;ind_int(ind_int>nbinfit)=nbinfit;              % edges intercept
    ind_grad(ind_grad<1)=1;ind_grad(ind_grad>nbinfit)=nbinfit;          % edges gradient
    
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
%---------Fit the straight line trend through the conditioned crossplot----------
% Fit y=mx straight line
if method ==1
    [p1,flag1]=lsqr(intercept_w,gradient_w);
    slope1=p1;g_int1=0;
    
    [p2,flag2]=lsqr(gradient_w,intercept_w);
    slope2=1/p2;g_int2=0;
end

% Fit y=mx +c straight line
if method==2
    p1=polyfit(intercept_w,gradient_w,1);% Solve minimizing  y-error
    %p=polyfit(gradient,intercept,1);
    slope1=p1(1);
    g_int1=p1(2);
    
    p2=polyfit(gradient_w,intercept_w,1);% Solve minimizing  x-error
    slope2=1/p2(1);
    g_int2=-p2(2)/p2(1);
end

slope_z=(slope1+slope2)/2;      % Average the two fits
grad0_z=(g_int1+g_int2)/2;


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
dist_trend=grad0_z*sind(chi_trend);
end

%% ########################### Function to return the peaks and troughs in a trace and their location as
% samples
function [sample_p sample_t trace_p trace_t]=peaks_and_troughs(trace)

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
%% ##################function to postprocess poorly conditioned data and replace it by median value of the result array
