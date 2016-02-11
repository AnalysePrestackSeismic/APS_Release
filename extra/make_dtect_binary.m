function [] = make_dtect_binary(input_dir, filename_pattern, n_files)
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

    cd(input_dir);
    
    start_time = 0;
    sample_rate = 4;

    for ii = 1:1:n_files
        load(sprintf('%s_%d.mat',filename_pattern,ii));
        if ii == 1
            n_samples = size(results_out{1,2},1);
        end
        fid_ilxl_f32 = fopen('ilxl_f32_tmp','w');
        fwrite(fid_ilxl_f32,results_out{5,2}{1},'uint32');
        fclose(fid_ilxl_f32);
        fid_ilxl_f32 = fopen('ilxl_f32_tmp','r');
        ilxl_f32_tmp = fread(fid_ilxl_f32,size(results_out{5,2}{1}),'*float32');
        fclose(fid_ilxl_f32);

        for kk = 1:1:size(results_out,1)-4
            fid_dtect_binary = fopen(sprintf('%s.bin',results_out{kk,1}),'a');
            if ii == 1
                fwrite(fid_dtect_binary,[start_time;sample_rate],'float32');
                fwrite(fid_dtect_binary,n_samples,'uint32');
            end
            fwrite(fid_dtect_binary,[ilxl_f32_tmp';results_out{kk,2}],'float32');
            fclose(fid_dtect_binary);
            if kk == size(results_out,1)-1
                fprintf('Written block %d of %d\n',ii,n_files);
            end
        end
    end   
end

