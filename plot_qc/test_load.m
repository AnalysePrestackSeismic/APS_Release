function test_load(job_meta_path,i_vol,i_block)
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
job_meta = load(job_meta_path);

[~, traces, ilxl_read] = node_segy_read(job_meta_path,i_vol,i_block);

results_out{1,1} = 'Meta data for output files';
results_out{1,2}{1,1} = ilxl_read;
results_out{1,2}{2,1} = uint32(zeros(size(traces,2)));
results_out{2,1} = 'test_load_write';
results_out{2,2} = traces;

node_segy_write(results_out,str2double(i_block),job_meta.s_rate/1000,job_meta.output_dir)
end