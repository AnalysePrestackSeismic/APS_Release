function [ output_args ] = near_far_calc(job_meta_path,i_block)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Load job meta information 
job_meta = load(job_meta_path);

% Load far volume
[~, far, ~] = ...
        node_segy_read(job_meta_path,'2',i_block);

[wb_idx] = water_bottom_picker(far,0);

% Load near volume
[~, near, ~] = ...
        node_segy_read(job_meta_path,'1',i_block);
% Time alignment and balancing


% Crossplot Near (N) versus Far Minus Near (FN) similar trend to intercept
% versus gradient
% Far minus Near times Far (FNXF) enhance class II
fn = far-near;
fnxf = fn.*far;

% Far minus Near times Near (FNXN) enhance class III
fnxn = fn.*near;

% Flatten to water bottom
[fn] = trace_flatten(fn,wb_idx);
[fnxf] = trace_flatten(fnxf,wb_idx);
[fnxn] = trace_flatten(fnxn,wb_idx);

% Transpose to make slice ordered
ns = size(fn,1);
fn = fn';
fnxf = fnxf';
fnxn = fnxn';

% Write as binary
fid_fn = strcat(job_meta.output_dir,'slice_fn_block',i_block,'.bin');
fid_fnxf = strcat(job_meta.output_dir,'slice_fnxf_block',i_block,'.bin');
fid_fnxn = strcat(job_meta.output_dir,'slice_fnxn_block',i_block,'.bin');
fwrite(fid_fn,ns,'float32');
fwrite(fid_fnxf,ns,'float32');
fwrite(fid_fnxn,ns,'float32');

fwrite(fid_fn,fn,'float32');
fclose(fid_fn);
fwrite(fid_fnxf,fnxf,'float32');
fclose(fid_fnxf);
fwrite(fid_fnxn,fnxn,'float32');
fclose(fid_fnxn);

end

