function [ freq_attr ]= avg_freq_attribute_calc (job_meta_path,st_vol,end_vol,blk,plot_pic,Z_diff,Qmax,power_cut_perc)
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
%% Function to Calculate Frequency Attributes and  a Q function based on Spectrums estimated in Wavelet Estimation
% This function calculates frequency based attributes from the estimated
% wavelets from Wavelet Estimation Algorithm
% Inputs:
% Outputs:
% Write to Disk:

%% ##############################################################################################################################################################
%-----Parameters------------
% Z_diff=1; % Rough desirabelDifference in the windows in z (time in s and depth in km);
% Qmax=250; % Maximum tolerated Q (This is not a hard threshold, its used for the Qestimation as a loose boundary

filt_len1=3;                                                                % Length of Smoothening Filter on Z axis in samples (used on wavelet matrix)
filt_len2=5;                                                                % Length of Smoothening Filter on frequency axis in samples (used on wavelet matrix)

Power_cut=power_cut_perc/100;                                                              % Fraction of maximum power (db)to find initial high and low boundaries on frequency scale

% This can be implemented by a smarter method using differential of the
% amplitude spectrum and looking for the point where it approaches zero
% using a tolerance



%% ###############################################################################################################################################################
%-----------Loading Stuff and creating Directories-------------
close all;
blk_struct=load('/data/TZA/segy/2013_pweza_chewa_presdm_pgs/gather_for_eage/digi_noQ/wavelets/all_wavelets_time.mat');
blk_cell=struct2cell(blk_struct);                           % erik: convert blk file to cell type
blk=cell2mat(blk_struct);                                   % erik: convert blk file to mat type
job_meta = load(job_meta_path);                             % Loading Job Meta file
wave_mat_path=strcat(job_meta.wav_directory,'smo_wavelet_algo_2_2_block_',num2str(blk)); % Path of Wavelet Mat file
wave = load(wave_mat_path);                                 %Loading Wavelet
n_vol = size(wave.all_wavelets_freq,2);                     % Number os angles or offsets in the wavelet file
[ns,n_win]=size((wave.all_wavelets_freq{1,2}));             % Number of Samples in the wavelets and number of wavelets
n_s=ns-1;
hn_s=n_s/2;                                                 % Half of Number of Samples
n_f= 1e6/job_meta.s_rate;
n_s_shift=job_meta.ns_win-job_meta.ns_overlap;
length_window=job_meta.s_rate*n_s_shift/1000;               % Z Length of the movement of moving window
n_diff = round(Z_diff*1000/length_window);                  % Rounder number of translation of the moving window to accomodate user supplied Z_diff

if mod(n_diff,2)==1
    n_diff=n_diff+1;                                        % Ensure odd number of samples between the compared windows so that the Q value can be assigned to the middle sample
end


Z_diff2=n_diff*length_window/1000;                          % Recalculated  Z_diff based of integer number of translations
n_vol= min(n_vol,abs(end_vol - st_vol));                    % Limiting number of volumes to user's specifications

if exist(strcat(job_meta.wav_directory,'Frequency_Attributes/'),'dir') == 0
    mkdir(strcat(job_meta.wav_directory,'Frequency_Attributes/'))
end

fig_dir = strcat(job_meta.wav_directory,'Frequency_Attributes/blk',num2str(blk),'_figs/');% Directory to Save Figures

if exist(fig_dir,'dir') == 0
    mkdir(fig_dir);
end
%% ###################################################################################################################################################################
%------------- Precondition wavelets/Spectrum---------------------

%-------Loading and QC of Wavelet Matrices------------------
w_freq=zeros(n_s,n_win,n_vol);                                              % Initialize matrix for wavelet frequency spectrum
w_time=zeros(n_s,n_win,n_vol);                                              % Initialize matrix for wavelet time series

for i_vol=st_vol:end_vol
    w_freq(:,:,i_vol) = wave.all_wavelets_freq{1,i_vol}(2:end,:);           % Filling matrix for wavelet frequency spectrum
    w_time(:,:,i_vol) = wave.all_wavelets_time{1,i_vol}(2:end,:);           % Filling matrix for wavelet time series
end

%-------------Check for NaNs-----------------------------------
%If freequency spectrum of an window has NaNs, copy the frequency spectrum
%of previous window to this window
% for i_vol=st_vol:end_vol
%     for k=2:n_win
%         if isnan(max(w_freq(:,k,i_vol))) > 0
%          w_freq(:,k,i_vol) = w_freq(:,k-1,i_vol) ;
%          w_time(:,k,i_vol) = w_time(:,k-1,i_vol) ;
%         end
%     end
% end

w_nan_flag=isnan(w_freq);                                                   % Find out the NaN is the Matrix (means there int a valid wavelet for these entries)
fold = sum (~w_nan_flag,3);                                                 % Count the number of valid wavelets along tertiar key (offset)
w_freq(isnan(w_freq))=0;                                                      % Replace the NaNs by zeros

% ----------COndition wavelets ---------------------------


w_freq_avg=bsxfun(@rdivide,sum(w_freq(:,:,st_vol:end_vol),3),fold);         % Average wavelet matrix across offset
w_freq_avg=w_freq_avg(1:hn_s,:);                                            % Retain only half part of it (since the spectrum is mirrored
w_freq_avg =20*log10(w_freq_avg);                                           % Converting wavelet matrix to db scale

%----------------Create axes------------------------------
freq_axis = (n_f)/2*linspace(0,1,job_meta.ns_win/2);                        % Create the frequency/wavenumer axis
z_axis= (job_meta.s_rate/1000)*(n_s/2:n_s_shift:n_s_shift*(n_win+1))';      % Create the z axis starting from middle of first window
close;

%-------------------------2D Filtering Wavelet Matrix-------------------
filt_design=ones(filt_len1,filt_len2);                                      % FIR filter to smoothen wavelet matrix


input_pad= w_freq_avg;                                                      % Initializing Input pad to the Wavelet Matrix
row_len=size(input_pad,2);                                                  % Number of rows in the wavelet matrix
col_len=size(input_pad,1);                                                  % Number of columns in the wavelet matrix

input_pad = [repmat(input_pad(:,1),1,(filt_len1+1)), input_pad, repmat(input_pad(:,row_len),1,(filt_len1+1))];      % Padding Rows by filter length
input_pad = [repmat(input_pad(1,:),(filt_len2+1),1); input_pad; repmat(input_pad(col_len,:),(filt_len2+1),1)];      % Padding Columns by filter length

input_pad_filt=conv2(input_pad,filt_design,'same');                                                                        % Smoothenig wavelet matrix
w_freq_avg_filt = input_pad_filt(((filt_len2+2):(filt_len2+2+col_len-1)),((filt_len1+2):(filt_len1+2+row_len-1)));  % Cropping the padded filtered matrix to original size
w_freq_avg_filt=w_freq_avg_filt/(filt_len1*filt_len2);                                                              % Dividing by number of samples being averaged over
clear input_pad row_len col_len;

%% ##########################################################################################################################################################################
% ---------Maximum Spectral Amplitude and Dominant Frequency Normalize all wavelets-----------
w_max_amp = max (w_freq_avg_filt, [],1);                                    % Maximum Spectral Amplitude
w_min_amp = min (w_freq_avg_filt, [],1);                                    % Maximum Spectral Amplitude
%scalar_norm= 1./w_max_amp;                                                  % Scalar for Normalization
w_dom_freq_ind = zeros(1,n_win);                                            % Initialize array for dominant frequency
w_dying_freq_ind = zeros(1,n_win);                                          % Initialize array for dying frequency
w_birth_freq_ind = zeros(1,n_win);                                          % Initialize array for birth frequency
w_freq_norm=zeros(n_s/2,n_win);                                             % Initializematrix dor normalized frequncy spectrum

for i=1:n_win
    
    w_dom_freq_ind(i)= find(w_freq_avg_filt(:,i)==w_max_amp(i));            % finding index of dominant frequency
    w_freq_norm(:,i)= (w_freq_avg_filt(:,i)-w_min_amp(i))/(w_max_amp(i)-w_min_amp(i));                  % Normalizing the spectrum between maximum and minimum
    
    w_dying_freq_ind(i)= find(w_freq_norm(:,i) > Power_cut,1,'last');       % finding index of dying frequency  (end limit of usable frequency spectrum)
    w_birth_freq_ind(i)= find(w_freq_norm(:,i) > Power_cut,1,'first');      % finding index of birth frequency (beginning limit of usable frequency spectrum)
    
end
w_dom_freq = freq_axis(w_dom_freq_ind);                                     % Calculate Dominant Frequency
w_dying_freq = freq_axis(w_dying_freq_ind);                                 % Calculate Dying Frequency
w_birth_freq = freq_axis(w_birth_freq_ind);                                 % Calculate Birth Frequency
bandwidth = w_dying_freq-w_birth_freq;                                      % Calculate Effective Bandwidth

clear scalar_norm i;
%% ####################################################################################################################################################################
%-----Calculating Slope of Decay---
S_decay_raw = 10*ones(1,n_win,'double');                                        % Initialize array for slope of decay of lee side of amplitude spectrum
C_decay_raw = zeros(1,n_win);
figure(6);
for  i = 1:n_win
    spec_amp2 = w_freq_avg_filt((w_dom_freq_ind(i)+0):w_dying_freq_ind(i),i);   % Cropping the ampitude spectrum in the right frequency range
    freq_axis2= freq_axis((w_dom_freq_ind(i)+0):w_dying_freq_ind(i));           % Crop relevant part of frequency axis
    p = polyfit (freq_axis2,spec_amp2',1);                                      % Find the slope of a linear fit between spectral amplitude and frequency
    S_decay_raw(i) = p(1);
    C_decay_raw(i) = p(2);
    
    %     plot(freq_axis2,spec_amp2,'b');
    %     hold on;
    %     plot(freq_axis2,S_decay_raw(i)*freq_axis2+C_decay_raw(i),'m');
    %     scatter ([w_dom_freq(i) w_dying_freq(i)],[30 30]);
    %     hold off;
    
end
hold off;
close all;


%% ############################################################################################################################################################################
% -----------Calculate Q---------------------------
close;
%Initialize the variables
Q =  zeros(1,n_win);                                                        % Initialize array for Attenuation from Conditioned Spectral Ratio
MUF= zeros(1,n_win);                                                        % Initialize array for Maxuimum Usable Frequency from Conditioned Spectral Ratio
Q_raw = zeros(1,n_win);                                                     % Initialize array for Attenuation from Raw Spectral Ratio
Dp = 10*ones(1,n_win,'double');                                             % Initialize array for slope from Spectral Ratio
Dp_raw = 10*ones(1,n_win,'double');                                         % Initialize array for slope from Raw Spectral Ratio
figure;
Dpmax=27.3*Z_diff2/Qmax;                                                    % Maximum Slope Tolerated (loose constraint)
S_decay = S_decay_raw;                                                      % Initialize array for slope of decay of lee side of amplitude spectrum from conditioned points
C_decay = C_decay_raw;

% Loop to calcuate Q in feed backward system
for i=(n_diff+1):n_win
    Spec_Ratio= w_freq_avg_filt((w_dom_freq_ind(i)):w_dying_freq_ind(i),i-n_diff)-w_freq_avg_filt((w_dom_freq_ind(i)):w_dying_freq_ind(i),i);   % Calculating Spectral Ratio of Amplitude
    freq_axis2= freq_axis((w_dom_freq_ind(i)):w_dying_freq_ind(i));                                                                             %Corresponding frequency axis
    
    Spec_Ratio_raw =w_freq_avg((w_dom_freq_ind(i)):w_dying_freq_ind(i),i-n_diff)-w_freq_avg((w_dom_freq_ind(i)):w_dying_freq_ind(i),i);         % Calculating Spectral Ratio of Amplitude
    Amp_Spec= w_freq_avg_filt((w_dom_freq_ind(i)):w_dying_freq_ind(i),i);                                                                       % Storing Amplitude Spectrum of the ith window
    %--------Preconditioning the spectrum and rarios-----------------
    Spec_Ratio2 = medfilt1(Spec_Ratio,3);                                   % Median Filter the Spectral Ratios along frequency axis
    Inputpad= [Spec_Ratio2(1);Spec_Ratio2(:)];                              % Pad the Spectral Ratio in the start
    diff_SpecRatio=diff(Inputpad);                                          % Differentiate the padded Spectral Ratio Array (for size matching)
    diff_SpecRatio(2)=diff_SpecRatio(1);                                    % Since the first entry will be zero duplicate the 2nd entry into the first (diffrence array)
    max_diff = max(diff_SpecRatio);                                         % Maximim Difference between consecutive frequency samples
    
    %--Rejection of sampes to ensure a stable positive slope
    k=1;                                                                    % Increment Index, incremented with in loop
    Spec_Ratio3(k) =Spec_Ratio2(1) ;                                        % Accept first sample from spectral ratio array
    freq_axis3(k) =freq_axis2(1) ;                                          % Accept first sample from frequency axis array
    Amp_Spec3 (k) = Amp_Spec(1);                                            % Accept first sample from spectral amplitude array
    
    %-Rejection Loop--
    for j=2:length(Spec_Ratio2)
        if diff_SpecRatio(j) >= (0.2*max_diff)                              % Reject if the the differential is lesser than 20% of max diffrential observed
            tol = (freq_axis2(j)-freq_axis2(1))*Dpmax+Spec_Ratio2(1);       % Minimum next Spectral Ratio value tolerated. If its greater it will be rejected
            if Spec_Ratio2(j) >tol;                                         % Reject sample if lesser than the claculated tolerance
                k=k+1;                                                      % Increment Index
                Spec_Ratio3(k) =Spec_Ratio2(j) ;                            % If all criteron satisfied accept the sample into a new array
                freq_axis3(k) =freq_axis2(j) ;                              % If all criteron satisfied accept corresponding frequency  into a new array
                Amp_Spec3 (k) = Amp_Spec(j);                                % If all criteron satisfied accept corresponding spectral amplitude  into a new array
            end
        end
    end
    
    %--------Plotting and line fiiting through the conditioned scatter points--------------
    %scatter(freq_axis3,Spec_Ratio3);%axis([10 100 -1 30]);                  % Scatter plot comment if you want
    %pause(1);
    if length(Spec_Ratio3)>1
        p = polyfit (freq_axis3,Spec_Ratio3,1);                             % Fitting a line through frequency and Spectral Ratios (conditioned)
        Dp(i-n_diff/2) = p(1);                                              % The slope term in the linear regresssion fit, This value is then is aligned to the sample in the middle of the compared Z windows
        p2 = polyfit (freq_axis3,Amp_Spec3,1);                              % Fitting a line through frequency and Spectral Amplitudes (conditioned)
        S_decay(i) = p2(1);                                                 % The slope term in the linear regresssion fit
        C_decay(i) = p2(2);                                                 % The intercept term in the linear regresssion fit
        %MUF(i-n_diff/2)=max(freq_axis3);                                    % The maximum usable frequency
        MUF(i)=max(freq_axis3);                                    % The maximum usable frequency
    else
        Dp(i-n_diff/2)=Dp(i-n_diff/2-1);                                    % If there's all points are rejected (The first point is accepted by default. Assign the previous Dp.
        S_decay(i)= S_decay(i-1);
        C_decay(i)= C_decay(i-1);
       % MUF(i-n_diff/2)=MUF(i-n_diff/2-1);
        MUF(i)=MUF(i);
    end
    Q(i-n_diff/2) = 27.3*Z_diff2/Dp(i-n_diff/2);                            %Calculate Q (effective over 1 sec data)
    
    
    %------------- Line fitting to raw unconditioned data points--------------
    p = polyfit (freq_axis2,Spec_Ratio_raw',1);                                 % Fitting a line through frequency and Spectral Ratios
    Dp_raw(i-n_diff/2) = p(1);
    Q_raw(i-n_diff/2) = 27.3*Z_diff2/Dp_raw(i-n_diff/2);                                %Calculate Q (effective over 1 sec data)
    clear Spec_Ratio Spec_Ratio2 Spec_Ratio3 freq_axis3 freq_axis2 k j p diff_SpecRatio max_diff Inputpad Amp_Spec3;
end
% Calculate Q, This value is thenis aligned to the sample in the middle of the compared Z windows
Dp(1:n_diff/2)= Dp(n_diff/2+1);                                             % Extrapolate for initial values
Q(1:n_diff/2)=Q(n_diff/2+1);                                                % Extrapolate for initial values
MUF(1:n_diff)=MUF(n_diff+1);                                            % Extrapolate for initial values
Q((n_win-n_diff/2):n_win)=Q(n_win-n_diff/2-1);                              % Extrapolate for final values
Dp((n_win-n_diff/2):n_win)=Dp(n_win-n_diff/2-1);                            % Extrapolate for final values
MUF((n_win-n_diff/2):n_win)=MUF(n_win-n_diff/2-1);                            % Extrapolate for final values

Dp_raw(1:n_diff/2)= Dp_raw(n_diff/2+1);                                     % Extrapolate for initial values
Q_raw(1:n_diff/2)=Q_raw(n_diff/2+1);                                        % Extrapolate for initial values
Q_raw((n_win-n_diff/2):n_win)=Q_raw(n_win-n_diff/2-1);                      % Extrapolate for final values
Dp_raw((n_win-n_diff/2):n_win)=Dp_raw(n_win-n_diff/2-1);                    % Extrapolate for final values

w_usable_bandwidth=MUF-w_birth_freq;                                          % Usable Bandwidth

close all;


%% #################################################################################################################################################################
%------------Plotting Results--------------------------------
close;
if plot_pic==1
    
    % Plot 1 : Spectrum
    h1=figure(1);
    subplot(2,3,1);plot(freq_axis,w_freq_avg(1:hn_s,:)); title('Frequency Spectrum');xlabel('Frequency');
    subplot(2,3,2);plot(freq_axis,w_freq_avg_filt(1:hn_s,:)); title('2D Smoothed Frequency Spectrum');xlabel('Frequency');
    subplot(2,3,3);plot(freq_axis,w_freq_norm(1:hn_s,:)); title('Normalized Frequency Spectrum');xlabel('Frequency');
    
    subplot(2,3,4);imagesc(freq_axis,z_axis,w_freq_avg(1:hn_s,:)');title('Frequency Spectrum');xlabel('Frequency');ylabel('Time / Depth');
    subplot(2,3,5);imagesc(freq_axis,z_axis,w_freq_avg_filt(1:hn_s,:)');title('2D Smoothed Frequency Spectrum');xlabel('Frequency');ylabel('Time / Depth');
    subplot(2,3,6);imagesc(freq_axis,z_axis,w_freq_norm(1:hn_s,:)');title('Normalized Frequency Spectrum');xlabel('Frequency');ylabel('Time / Depth');
    
    saveas(h1,strcat(fig_dir,'Spectrum_blk',num2str(blk)),'png');
    %close all;
    
    % Plot 2: Results
    h2=figure (3);
    %----------------------------------------------------------------------------
    subplot(1,5,1);plot(w_max_amp,z_axis,'-r','LineWidth',2); set(gca,'YDir','rev');
    title('Maximum Spectral Amplitude');ylabel('depth/time');xlabel('Amplitude');ylim([min(z_axis) max(z_axis)]);
    
    %----------------------------------------------------------------------------
    subplot(1,5,2);plot(w_dom_freq,z_axis,'-b','LineWidth',2);set(gca,'YDir','rev');
    title('Frequencies');ylabel('depth/time');xlabel('Frequency');ylim([min(z_axis) max(z_axis)]);xlim([0 100]);
    hold all;
    
    subplot(1,5,2);plot(w_birth_freq,z_axis,'-m','LineWidth',2);set(gca,'YDir','rev');
    title('Frequencies');ylabel('depth/time');xlabel('Frequency');ylim([min(z_axis) max(z_axis)]);xlim([0 100]);
        
    hold all;
    subplot(1,5,2);plot(w_dying_freq,z_axis,'-g','LineWidth',2);set(gca,'YDir','rev');
    title('Frequencies');ylabel('depth/time');xlabel('Frequency');ylim([min(z_axis) max(z_axis)]);xlim([0 100]);
    
    hold all;
    subplot(1,5,2);plot(MUF,z_axis,'-r','LineWidth',2);set(gca,'YDir','rev');
    title('Frequencies');ylabel('depth/time');xlabel('Frequency');ylim([min(z_axis) max(z_axis)]);xlim([0 100]);
    %---------------------------------------------------------------------------
    
    subplot(1,5,3);plot(bandwidth,z_axis,'-c','LineWidth',2);set(gca,'YDir','rev');
    title('Effective Bandwidth');ylabel('depth/time');xlabel('Frequency');ylim([min(z_axis) max(z_axis)]);xlim([10 100]);
    
    hold all;
    
    subplot(1,5,3);plot(w_usable_bandwidth,z_axis,'--m','LineWidth',2);set(gca,'YDir','rev');
    title('Effective Bandwidth');ylabel('depth/time');xlabel('Frequency');ylim([min(z_axis) max(z_axis)]);xlim([10 100]);
    %---------------------------------------------------------------------------
    
    subplot(1,5,4);plot(S_decay_raw,z_axis,'-b','LineWidth',2);set(gca,'YDir','rev');
    title('Slope of Decay');ylabel('depth/time');xlabel('Attenuation (Qeff)');ylim([min(z_axis) max(z_axis)]);xlim([-1 0]);
    hold all;
    
    subplot(1,5,4);plot(S_decay,z_axis,'-r','LineWidth',2);set(gca,'YDir','rev');
    title('Slope of Decay');ylabel('depth/time');xlabel('Attenuation (Qeff)');ylim([min(z_axis) max(z_axis)]);xlim([-1 0]);
    %---------------------------------------------------------------------------
    
    subplot(1,5,5);plot(Q,z_axis,'-m','LineWidth',2);set(gca,'YDir','rev');
    title('Effective Attenuation');ylabel('depth/time');xlabel('Attenuation (Qeff)');xlim([0 300]);ylim([min(z_axis) max(z_axis)]);
    hold all;
    
    subplot(1,5,5);plot(Q_raw,z_axis,'-g','LineWidth',2);set(gca,'YDir','rev');
    title('Effective Attenuation Raw');ylabel('depth/time');xlabel('Attenuation (Qeff)');xlim([0 300]);ylim([min(z_axis) max(z_axis)]);
    hold off;
    
    saveas(h2,strcat(fig_dir,'Output_blk',num2str(blk)),'png');
    %close all;
    
    % PlotQC
    h3=figure(4);
    l=1;
    for win = floor(n_win/5):floor(n_win/5):n_win;
        
        subplot(1,5,l);
        a_cond =S_decay(win)*freq_axis+C_decay(win);
        a_raw = S_decay_raw(win)*freq_axis+C_decay_raw(win);
        
        plot(freq_axis,w_freq_avg(1:hn_s,win),'g'); title('Frequency Spectrum');xlabel('Freuency');xlim([0 100]);ylim([0 100])
        hold on;
        plot(freq_axis,w_freq_avg_filt(1:hn_s,win),':m'); title('Frequency Spectrum');xlabel('Freuency');xlim([0 100]);ylim([0 100])
        
        plot(freq_axis,a_cond,'b'); title('Frequency Spectrum');xlabel('Frequency');
        
        plot(freq_axis,a_raw,'r'); title('Frequency Spectrum');xlabel('Frequency');
        hold off;
        l=l+1;
    end
    
    saveas(h3,strcat(fig_dir,'QC_Slope_blk',num2str(blk)),'png');
    %close all;
end
clear h1 h2 h3 l i loc length_window;
%end

%% ######################################################################################################################################################################
% Writting Output to Disk

loc= job_meta.block_keys(blk,:);
inl_central= floor((loc(1)+loc(2))/2);
xl_central = floor((loc(3)+loc(4))/2);
inl_ary = inl_central*ones(length(z_axis),1);
xl_ary = xl_central*ones(length(z_axis),1);
name = blk*ones(length(z_axis),1);

freq_attr = struct ('inl_central',inl_central,'xl_central',xl_central,'z_axis',z_axis,'max_spec_amp',w_max_amp,'peak_freq',w_dom_freq,'start_freq',w_birth_freq,'end_freq',w_dying_freq,'effective_end_freq',MUF,'s_decay',S_decay,'q_conditioned',Q,'q_raw',Q_raw);
freq_attr_Mat =[name z_axis inl_ary xl_ary w_max_amp' w_dom_freq' w_dom_freq' w_birth_freq' w_dying_freq' S_decay' Q' Q_raw'];

file_str = strcat(job_meta.wav_directory,'Frequency_Attributes/','Frequency_Attribute',num2str(blk),'.txt');

fid_wav = fopen(file_str ,'w+');
header = strcat('Frequency Attributes for Block Number: ', num2str(blk),sprintf('\n%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n','--------------------------------','Z Axis (ms/m)','Max_Spec_Amp','Peak_Freq','start_Freq','end_freq','S_decay','Q_conditioned','Q_raw','-------------------------------'));

fprintf(fid_wav, '%s\n',header);
fprintf(fid_wav, '%s%g\n','Central Inline : ',inl_central);
fprintf(fid_wav, '%s%g\n','Central Xline : ',xl_central);


fclose(fid_wav);
dlmwrite(file_str,freq_attr_Mat,'delimiter', '\t','precision', '%4.0f' ,'-append');

end

