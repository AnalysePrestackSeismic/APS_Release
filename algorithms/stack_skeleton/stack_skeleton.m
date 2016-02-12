function [ skeleton ] = stack_skeleton( inputsegy,varargin )
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
%STACK_SKELETON Generate a binary mask from a stack with ones at the
%peaks/troughs of all spatially continuous events
%
%   inputsegy:      segy file of input stack
%   optional arguments:
%   'ilbyte='           inline byte location (defaults to 189)
%   'xlbyte='           xline byte location (defaults to 193)
%   'threshold='        threshold for picking peaks, relative to overall
%                       rms amplitude (default 1.5)
%   'segment='          minimum segment size in pixels (default 100)
%   'gap='              gap between peaks, in samples (defaults to 5)

% Defaults for arguments =================================================
ilbyte = 189;
xlbyte = 193;
threshold = 1.5;
gap = 10;
segment = 100;
%=========================================================================
%
for kv = 1:length(varargin)
    varknown = false;
    if strfind(varargin{kv},'ilbyte=')
        vartmp = deblank(regexp(varargin{kv},'=','split'));
        ilbyte = str2double(vartmp(2));
        varknown = true;
    end
    if strfind(varargin{kv},'xlbyte=')
        vartmp = deblank(regexp(varargin{kv},'=','split'));
        xlbyte = str2double(vartmp(2));
        varknown = true;
    end
    if strfind(varargin{kv},'threshold=')
        vartmp = deblank(regexp(varargin{kv},'=','split'));
        threshold = str2double(vartmp(2));
        varknown = true;
    end
    if strfind(varargin{kv},'gap=')
        vartmp = deblank(regexp(varargin{kv},'=','split'));
        gap = str2double(vartmp(2));
        varknown = true;
    end
     if strfind(varargin{kv},'segment=')
        vartmp = deblank(regexp(varargin{kv},'=','split'));
        segment = str2double(vartmp(2));
        varknown = true;
    end
   
    % last test to see if there was a command line variable that was not
    % expected in which case throw warning and crash
    if varknown == false;
        error('unknown input variable; the variable %s is not recognised and the function has quit\n', varargin{kv});
    end
end

% read the stack traces

[stack_str stack_ilxl stack_traces]=segy_to_mat(ilbyte,xlbyte,inputsegy);

% balance amplitudes

stack_traces = time_balence(stack_traces);
% stack_traces = abs(stack_traces);

% stack_smth = gaussian_2dsmth(grid_stack(:,:,5),11,5); % replace with 3d smoothing

% get the RMS amplitude

rms = sqrt(mean(nonzeros(stack_traces).^2));
threshold = threshold * rms;

% reduce each trace to just the peak/trough values

stack_peaks = zeros(size(stack_traces));
stack_troughs = zeros(size(stack_traces));

for tt=1:stack_str.n_traces
    [stack_peaks(:,tt),~] = find_peaks(stack_traces(:,tt),threshold,gap);
end

for tt=1:stack_str.n_traces
    [stack_troughs(:,tt),~] = find_peaks(-stack_traces(:,tt),threshold,gap);
end

% put it into a grid

finl = min(stack_ilxl(:,1));
linl = max(stack_ilxl(:,1));
fxl = min(stack_ilxl(:,2));
lxl = max(stack_ilxl(:,2));

ilinc = min(nonzeros(diff(sort(stack_ilxl(:,1)))));
xlinc = min(nonzeros(diff(sort(stack_ilxl(:,2)))));

[ilxl_idx,ninl,nxl] = grid_index(stack_ilxl(:,1),stack_ilxl(:,2),finl,linl,ilinc,fxl,lxl,xlinc);

grid_peaks = zeros(stack_str.n_samples,ninl*nxl);
grid_peaks(:,ilxl_idx) = stack_peaks(:,:);
grid_peaks = reshape(grid_peaks,stack_str.n_samples,nxl,ninl);

grid_troughs = zeros(stack_str.n_samples,ninl*nxl);
grid_troughs(:,ilxl_idx) = stack_troughs(:,:);
grid_troughs = reshape(grid_troughs,stack_str.n_samples,nxl,ninl);



grid_peaks_con = bwareaopen(grid_peaks,segment,26);
grid_troughs_con = bwareaopen(grid_troughs,segment,26);

skeleton =grid_peaks_con + grid_troughs_con;
skeleton(skeleton~=0)=1;

figure; imagesc(squeeze(skeleton(:,:,5))); colormap('gray'); caxis([-1 1]);

% still have two big issues:
% 1. there are gaps along some events
% 2. we have events too close together because we're picking peaks and
% troughs


end

