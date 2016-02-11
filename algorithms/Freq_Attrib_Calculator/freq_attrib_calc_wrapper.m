function [] = freq_attrib_calc_wrapper(job_meta_path,Z_diff,Qmax,power_cut_perc)
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
% -------------------------------------------------------------------------
% Function to calculate frequency attributes based on spectral nature of
% wavelets and therir variability in 3D. Currently this code doenst tray to
% estimate QVO

% INPUTS:
%   job_meta_path : the job mets path (note this code must be run after wavelet estimation is done)
%   Z_diff: The difference between two windows for doing spectral ratio
%   Q max : The maximum allowable value for Q (300 is a good guess)
%   power_cut_perc : The percentage of max power to find the initial
%   usuable bandwidth of the spectrum

% WRITE TO DISK:
%   
% -------------------------------------------------------------------------
plot_pic=0;
%--------------------------------------------------------------------------


job_meta = load(job_meta_path);                                                         % Load job meta information
s_rate = job_meta.s_rate/1000;                                                          % Sampling rate in ms
%z_start =s_rate* round(job_meta.wb_z_avg(job_meta.liveblocks)/s_rate);
n_s_shift=job_meta.ns_win-job_meta.ns_overlap;

ns_eff = max(cell2mat(job_meta.n_samples))/n_s_shift;                                   % Effective Number of samples

if job_meta.is_gather == 0                                                              % For the case of using angle stacks
    i_vol_max = job_meta.nvols;                                                         % Number of angle stacks
else                                                                                    % For the case of using angle gathers
    i_vol_max = ((job_meta.tkey_max -  job_meta.tkey_min )/ job_meta.tkey_inc) +  1;    % Number of angles in angle gather
end
output_dir=strcat(job_meta.wav_directory,date,'Frequency_Attribute_Volumes','_',num2str(Z_diff),'_',num2str(Qmax),'_',num2str(power_cut_perc),'/');

if exist(strcat(output_dir),'dir') == 0
    mkdir(output_dir);% Make a directory for storing the volumes
end



s{1}=strcat(output_dir,'max_spec_amp.txt');
s{2}=strcat(output_dir,'peak_freq.txt');
s{3}=strcat(output_dir,'start_freq.txt');
s{4}=strcat(output_dir,'end_freq.txt');
s{5}=strcat(output_dir,'s_decay.txt');
s{6}=strcat(output_dir,'q_conditioned.txt');
s{7}=strcat(output_dir,'q_raw.txt');
s{8}=strcat(output_dir,'effective_end_freq');



% Open files
fid_wave=zeros(1,8);
for i=1:8
    fid_wave(i) = fopen( s{i} ,'w+');
    fclose(fid_wave(i));
end
%  fid_wav2 = fopen( strcat(job_meta.wav_directory,'Frequency_Attribute_Volumes/peak_freq.txt') ,'w+');
%  fid_wav3 = fopen( strcat(job_meta.wav_directory,'Frequency_Attribute_Volumes/start_freq.txt'),'w+');
%  fid_wav4 = fopen( strcat(job_meta.wav_directory,'Frequency_Attribute_Volumes/end_freq.txt'),'w+');
%  fid_wav5 = fopen( strcat(job_meta.wav_directory,'Frequency_Attribute_Volumes/s_decay.txt'),'w+');
%  fid_wav6 = fopen( strcat(job_meta.wav_directory,'Frequency_Attribute_Volumes/q_conditioned.txt'),'w+');
%  fid_wav7 = fopen( strcat(job_meta.wav_directory,'Frequency_Attribute_Volumes/q_raw.txt'),'w+');

%fprintf(fid_wav1, '%s\n','test');fprintf(fid_wav2, '%s\n','test');
%fclose(fid_wav1); fclose(fid_wav2);fclose(fid_wav3); fclose(fid_wav4);fclose(fid_wav5);fclose(fid_wav6); fclose(fid_wav7);

%fprintf(fid_wav1,num2str(z_start),

%% Wavelet Database Management for all bocks-- preconditioning for averaging

loopfin = size(job_meta.liveblocks,1);                                      % Number of Live Blocks
lpi = 1;                                                                    % Loop Index
%count = 1;                                                                 % Counter for blocks that have a corresponding wavelet file
z_start = 0;
Z_s_rate = 0;
n_Zsamples = 0;
while lpi <= loopfin
    i_block = job_meta.liveblocks_il_sorted(lpi);                           % Reference the block numbers for live blocks (from innline sorted array)
    %fprintf('keys %g \n', job_meta.block_keys(i_block,:));
    
    
    freq_attr = avg_freq_attribute_calc (job_meta_path,1,i_vol_max,i_block,plot_pic,Z_diff,Qmax,power_cut_perc);   %Call the freequency attribute calculator for assessed block
    %--------------------Write to Ascii file----------------------------------------
    if(lpi==1)
        Z_start = freq_attr.z_axis(1);                                      % Beginning of Z axis
        Z_s_rate =freq_attr.z_axis(2) -freq_attr.z_axis(1);                 % Z sample rate (m or ms)
        n_Zsamples =60; %length(freq_attr.z_axis);                               % Number of Z samples
        header = [Z_start Z_s_rate n_Zsamples];
        for k=1:7
            dlmwrite(s{k},header,'delimiter',' ','precision', '%4.0f','newline','unix' ,'-append');% Write the start z , sampling rate and number of samples into the files
        end
        
    end
    
    fprintf('Assessing Block %g: il : %g , xl : %g  \n',i_block,freq_attr.inl_central,freq_attr.xl_central);                                % Display the current block being assessed   job_meta.block_keys(i_block,:)
    
    for k=1:8
        switch k
            case 1
                w_array = [freq_attr.inl_central freq_attr.xl_central freq_attr.max_spec_amp];  % Build the trace preceded by position (il xl)
            case 2
                w_array = [freq_attr.inl_central freq_attr.xl_central freq_attr.peak_freq];     % Build the trace preceded by position (il xl)
            case 3
                w_array = [freq_attr.inl_central freq_attr.xl_central freq_attr.start_freq];    % Build the trace preceded by position (il xl)
            case 4
                w_array = [freq_attr.inl_central freq_attr.xl_central freq_attr.end_freq];      % Build the trace preceded by position (il xl)
            case 5
                w_array = [freq_attr.inl_central freq_attr.xl_central freq_attr.s_decay];       % Build the trace preceded by position (il xl)
            case 6
                w_array = [freq_attr.inl_central freq_attr.xl_central freq_attr.q_conditioned]; % Build the trace preceded by position (il xl)
            case 7
                w_array = [freq_attr.inl_central freq_attr.xl_central freq_attr.q_raw];         % Build the trace preceded by position (il xl)
            case 8
                w_array = [freq_attr.inl_central freq_attr.xl_central freq_attr.effective_end_freq];         % Build the trace preceded by position (il xl)
        end
        
        w_array2=zeros(1,n_Zsamples+2);
        w_array2(1:length(w_array))=w_array;
        dlmwrite(s{k},w_array2,'delimiter',' ','precision', '%5.0f','newline','unix','-append'); %Write the trace in the file
    end
    
    %------------------------------end writing------------------------------
    
    lpi = lpi + 1;                                          % Increment loop index
end



end
