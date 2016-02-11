function [] = ava_analysis_wrapper(job_meta_path)
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
% AVA ANALYSIS WRAPPER

% INPUT: job meta path : this will access the last cpa directory writen in
% the

% OUTPUT:

% WRITES TO DISK :
% -------------------------------------------------------------------------
%%

job_meta = load(job_meta_path);                                                         % Load job meta information
inline_sort_blocks(job_meta_path)                                                       % Sort block index in inline sorted fashion

output_dir=strcat(job_meta.cpa_directory,'cpa_joined_result/');


if exist(strcat(output_dir),'dir') == 0
    mkdir(output_dir);% Make a directory for storing the volumes
end


%%

loopfin = size(job_meta.liveblocks,1);                                      % Number of Live Blocks
lpi = 1;                                                                    % Lopp index 1
lpi2= 1;                                                                   % Loop Index (this is more like a flag)
while lpi <= loopfin
    i_block = job_meta.liveblocks_il_sorted(lpi);                           % Reference the block numbers for live blocks (from innline sorted array)
    
    ava_file_path = strcat(job_meta.cpa_directory,'crossplot_analysis_block',num2str(i_block),'.mat');
    
    if exist(ava_file_path,'file')~=0
        
        
        results_AVA=load(ava_file_path);   %Call the frequency attribute calculator for assessed block
        %--------------------Write to Ascii file----------------------------------------
        
        % If this is the first block accessed
        if(lpi2==1)
            %----initilize the folders
            s=results_AVA.header;
            n_items=size(results_AVA.data,1);   % Find the number of items to be written out
          
            Z_start = results_AVA.zaxis(1);                                      % Beginning of Z axis
            Z_s_rate =results_AVA.zaxis(2) -results_AVA.zaxis(1);                 % Z sample rate (m or ms)
            n_Zsamples =length(results_AVA.zaxis);                               % Number of Z samples
            header = [Z_start Z_s_rate n_Zsamples];
            for k=1:n_items
                dlmwrite(strcat(output_dir,s{k},'.txt'),header,'delimiter',' ','precision', '%4.0f','newline','unix');                  % Background:  Write the start z , sampling rate and number of samples into the files. Overwrite previous content
                dlmwrite(strcat(output_dir,s{k},'_extremas_','.txt'),header,'delimiter',' ','precision', '%4.0f','newline','unix');     % Peaks and troughs: Write the start z , sampling rate and number of samples into the files. Overwrite previous content
                dlmwrite(strcat(output_dir,s{k},'_peaks_','.txt'),header,'delimiter',' ','precision', '%4.0f','newline','unix');        % Peaks: Write the start z , sampling rate and number of samples into the files. Overwrite previous content
                dlmwrite(strcat(output_dir,s{k},'_troughs_','.txt'),header,'delimiter',' ','precision', '%4.0f','newline','unix');      % Troughs: Write the start z , sampling rate and number of samples into the files. Overwrite previous content
                
            end
            lpi2=lpi2-1;
        end
        
        % For all other blocks accessed
        for p=1:n_items
            w_array = [results_AVA.inl_central results_AVA.xl_central results_AVA.data(p,:)];
            w_array_pts=[results_AVA.inl_central results_AVA.xl_central results_AVA.data_pts(p,:)];
            w_array_p=[results_AVA.inl_central results_AVA.xl_central results_AVA.data_p(p,:)];
            w_array_t=[results_AVA.inl_central results_AVA.xl_central results_AVA.data_t(p,:)];
            
            
            %         w_array2=zeros(1,n_Zsamples+2);
            %
            %         w_array2(3:(length(w_array)+2))=w_array;
            dlmwrite(strcat(output_dir,s{p},'.txt'),w_array,'delimiter',' ','precision', '%4.0f','newline','unix','-append'); % Write the trace in the file. Append previous content
            dlmwrite(strcat(output_dir,s{p},'_extremas_','.txt'),w_array_pts,'delimiter',' ','precision', '%4.0f','newline','unix','-append');  % Peaks and Troughs together
            dlmwrite(strcat(output_dir,s{p},'_peaks_','.txt'),w_array_p,'delimiter',' ','precision', '%4.0f','newline','unix','-append');   % 
            dlmwrite(strcat(output_dir,s{p},'_troughs_','.txt'),w_array_t,'delimiter',' ','precision', '%4.0f','newline','unix','-append');
        end
        fprintf('Assessing Block %g: il : %g , xl : %g ||| Percentage complete %g \n',i_block,results_AVA.inl_central,results_AVA.xl_central, floor(100*(lpi/loopfin)));                                % Display the current block being assessed   job_meta.block_keys(i_block,:)
        %------------------------------end writing------------------------------
    end
    
    lpi = lpi + 1;                                          % Increment loop index
end
end
