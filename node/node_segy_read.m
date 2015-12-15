function [seismic, traces, ilxl_read, offset_read] = node_segy_read(job_meta_path,vol_index,i_block,varargin)
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
% node_segy_read: function to read traces from a specific block with a
% scanned segy volume
%   Arguments:
%       job_meta_path = path to job_meta .mat file
%       vol_index = integer indicating volume to load
%       i_block = integer indicating which block to load from volume
%       vargin: Optional extra flags for zmin and zmax to limit the
%               trace length 
%   Outputs:
%       seismic = structure containing seismic header information
%       traces = seismic traces as matrix (rows samples; columns traces)
%       ilxl_read = inlines and crosslines that have been read
%       offset_read = offsets that have been read
%
%   Writes to Disk:
%       nothing

%%
i_block = str2double(i_block);
job_meta = load(job_meta_path);

ztrunc = 0;
if ~isempty(varargin)
    if length(varargin) == 2
        zmin = varargin{1};
        zmax = varargin{2};
        ztrunc = 1;
    else
        error('must specify zmin and zmax on command line if using z limits');
    end
end

% calculate fold of gather - this should be added to job_meta
if job_meta.is_gather == 1
    fold = ((job_meta.tkey_max -  job_meta.tkey_min )/ job_meta.tkey_inc) +  1;
end 

% Error checking for standalone mode.
if i_block > job_meta.n_blocks
    seismic = NaN;
    traces = NaN;
    ilxl_read = NaN;
    offset_read = NaN;
    fprintf('Job only has %d blocks...\n',job_meta.n_blocks)
    return
end

% Lookup traces and find byte location in respective segy file
vol_keys = segy_index_byte_finder(job_meta_path,job_meta.block_keys(i_block,:),vol_index);
vol_index = str2double(vol_index);
vol_name = job_meta.volumes{vol_index};
    
% find the names of the blocks for the given volume
ii = 1;

if job_meta.is_gather == 0
    loop_index = size(job_meta.files,1);
else
    loop_index = size(job_meta.files,2);
end

for i_files = 1:1:loop_index
    if strfind(job_meta.files{i_files},vol_name) == 1
        blocks{ii,1} = job_meta.files{i_files};
        ii = ii + 1;
    end
end    
   
%% loop over the blocks
for ii_block = 1:1:size(blocks,1)
    if size(vol_keys{ii_block},1) > 2
        vol_keys{ii_block}(end+1,3) = 0;
        seismic = segy_read_binary(strcat(job_meta.paths{1},blocks{ii_block}));  

        bytes_per_sample = job_meta.bytes_per_sample{vol_index};
        trc_head = 240;
        trc_length = seismic.n_samples*bytes_per_sample; % could add this to seismic structure
        
        s_key = 1;
        is_key = 1;
        i_counter = 1;
        is_counter = 1;
        i_trace = 1;
        % find the number of rows in the bytes lookuptable
        while_is_gather = size(vol_keys{ii_block},1);
        % loop round all the row in the lookup table
        while s_key < while_is_gather
            if job_meta.is_gather == 0
                if vol_keys{ii_block}(i_counter+1,3) - vol_keys{ii_block}(i_counter,3) == trc_length+trc_head
                    i_counter = i_counter + 1;
                else  % perform continuous read
                    n_traces_to_read = i_counter-s_key+1;
                    % Open the seismic segy file
                    [traces{ii_block}(:,s_key:s_key+n_traces_to_read-1),...
                        ilxl_read{ii_block}(s_key:s_key+n_traces_to_read-1,:),...
                        offset_read{ii_block}(s_key:s_key+n_traces_to_read-1,:)] ...
                        = read_traces_segy(seismic,vol_keys{ii_block}(s_key,3)-trc_head,n_traces_to_read);
                    s_key = s_key + n_traces_to_read;
                    i_counter = s_key;
                end
            elseif job_meta.is_gather == 1
                % test to see if the gather is the same size as the
                % next one - this should see if you can do a continious
                % read
%                 if vol_keys{ii_block}(i_counter+1,3) - vol_keys{ii_block}(i_counter,3) == vol_keys{ii_block}(i_counter,4)*(trc_length+trc_head)
%                     i_counter = i_counter + 1;
%                     %is_counter = is_counter + 1;
%                 else  % perform continuous read
%                     n_traces_to_read = i_counter-s_key;
%                     %is_traces_to_read = is_counter-is_key+1;                    
%                     is_traces_to_read = (n_traces_to_read)*vol_keys{ii_block}(i_counter-1,4);
%                     % read to a temp array and then put into the correct
%                     % location in the output array
%                     
%                     % Open the seismic segy file
%                     [traces{ii_block}(:,is_key:is_key+is_traces_to_read-1), ilxl_read{ii_block}(is_key:is_key+is_traces_to_read-1,:), offset_read{ii_block}(is_key:is_key+is_traces_to_read-1,:)] ...
%                         = read_traces_segy(seismic,vol_keys{ii_block}(s_key,3)-trc_head,is_traces_to_read);    
%                     
%                     
%                     %alloffpres = size(offset_read{28},1)/fold
%                     
%                     s_key = s_key + n_traces_to_read;
%                     %i_counter = i_counter + 1;
%                     is_key = is_key + is_traces_to_read;
%                     %is_counter = is_counter + 1;
%                     end 
                
                %set the first gather to its own fold
                if i_counter == 1
                    prevfoldck = vol_keys{ii_block}(i_counter,4);
                end 
                    
                if vol_keys{ii_block}(i_counter+1,3) - vol_keys{ii_block}(i_counter,3) == vol_keys{ii_block}(i_counter,4)*(trc_length+trc_head) && vol_keys{ii_block}(i_counter,4) == prevfoldck
                    prevfoldck = vol_keys{ii_block}(i_counter,4);
                    i_counter = i_counter + 1;
                    is_counter = is_counter + 1;
                    %
                else  % perform continuous read
                    %set the first gather to its own fold
                    if i_counter == 1
                        prevfoldck = vol_keys{ii_block}(i_counter,4);
                    else
                        prevfoldck = vol_keys{ii_block}(i_counter-1,4);
                    end
                    %
                    if vol_keys{ii_block}(i_counter,4) == prevfoldck
                        n_traces_to_read = i_counter-s_key+1;
                    else
                        if i_counter == s_key
                            i_counter = i_counter + 1;
                            %n_traces_to_read = 1;
                        end
                        n_traces_to_read = i_counter-s_key;
                        
                    end
                    %is_traces_to_read = is_counter-is_key+1;
                    if i_counter == 1
                        is_traces_to_read = (n_traces_to_read)*vol_keys{ii_block}(i_counter,4);
                        curtestfold = vol_keys{ii_block}(i_counter,4);
                    else
                        is_traces_to_read = (n_traces_to_read)*vol_keys{ii_block}(i_counter-1,4);
                        curtestfold = vol_keys{ii_block}(i_counter-1,4);
                    end
                    % read to a temp array and then put into the correct
                    % location in the output array
                    
                    % Open the seismic segy file
%                     [traces{ii_block}(:,is_key:is_key+is_traces_to_read-1), ilxl_read{ii_block}(is_key:is_key+is_traces_to_read-1,:), offset_read{ii_block}(is_key:is_key+is_traces_to_read-1,:)] ...
%                         = read_traces_segy(seismic,vol_keys{ii_block}(s_key,3)-trc_head,is_traces_to_read);    
                    %[traces{ii_block}(:,is_key:is_key+is_traces_to_read-1), ilxl_read{ii_block}(is_key:is_key+is_traces_to_read-1,:), offset_read{ii_block}(is_key:is_key+is_traces_to_read-1,:)] ...
                    %    = read_traces_segy(seismic,vol_keys{ii_block}(s_key,3)-trc_head,is_traces_to_read);    
                    
%                     for cj = 1:10
%                         fprintf('%-20d%-20d%-20d%-20d\n',vol_keys{ii_block}(cj,1),vol_keys{ii_block}(cj,2),vol_keys{ii_block}(cj,3),vol_keys{ii_block}(cj,4))
%                     end
                    %alloffpres = size(offset_read{28},1)/fold ;% fprintf('%-20d\n',vol_keys{ii_block}(s_key,3))
                    
                    if fold == curtestfold
                        [traces{ii_block}(:,is_key:is_key+is_traces_to_read-1), ilxl_read{ii_block}(is_key:is_key+is_traces_to_read-1,:), offset_read{ii_block}(is_key:is_key+is_traces_to_read-1,:)] ...
                        = read_traces_segy(seismic,vol_keys{ii_block}(s_key,3)-trc_head,is_traces_to_read);   
                    else
                        [tmptracesblk, tmp_ilxl_read, tmp_offset_read] = read_traces_segy(seismic,vol_keys{ii_block}(s_key,3)-trc_head,is_traces_to_read);  
                        
                        if sum(diff(tmp_offset_read)) ~= 0
                            offset_inc=mode(diff(tmp_offset_read));
                        %if offset_inc==0
                        else
                            offset_inc= tmp_offset_read(1);
                        end
                        
                        %filltraces = zeros(size(tmptracesblk,1),(is_counter-is_key+1)*fold,'single');
                        
                        %traces{ii_block}(:,is_key:is_key+((is_counter-is_key+1)*fold)-1) = zeros(size(tmptracesblk,1),(is_counter-is_key+1)*fold,'single');
                        %ilxl_read{ii_block}(is_key:is_key+((is_counter-is_key+1)*fold)-1,:) = zeros((is_counter-is_key+1)*fold,size(tmp_ilxl_read,2),'int32');
                        
                        traces{ii_block}(:,is_key:is_key+((n_traces_to_read)*fold)-1) = zeros(size(tmptracesblk,1),(n_traces_to_read)*fold,'single');
                        ilxl_read{ii_block}(is_key:is_key+((n_traces_to_read)*fold)-1,:) = zeros((n_traces_to_read)*fold,size(tmp_ilxl_read,2),'int32');                        
                        
                        
                        % do unique on the tmp_ilxl_gather to get all the
                        % values
                        
                        %offset_read{ii_block}(is_key:is_key+((is_counter-is_key+1)*fold)-1,:) = zeros((is_counter-is_key+1)*fold,size(tmp_offset_read,2));
                        offset_read{ii_block}(is_key:is_key+((n_traces_to_read)*fold)-1,:) = repmat((int32(job_meta.tkey_min):int32(job_meta.tkey_inc):int32(job_meta.tkey_max))',(n_traces_to_read),1);
                        % ((job_meta.tkey_max -  job_meta.tkey_min )/ job_meta.tkey_inc) 
                        % the zeros
                        
                        
                        lastoffval = 99999999;
                        locval = -1;
                        finalloc = 1;
                        lastilxlval = -9876543;
                        listi = 1;
                        for cji = 1:1:is_traces_to_read;
                            if tmp_offset_read(cji) <= lastoffval && sum(tmp_ilxl_read(cji,:)) ~= lastilxlval;
                                locval = locval + 1;
                                tmpilxllist(listi,:) = tmp_ilxl_read(cji,:);
                                listi = listi + 1;
                            end
                            
%                             if isfield(job_meta,irreg_gath)
%                                                                 if job_meta.irreg_gath  ~= 1
%                                 finalloc = (locval * fold) + tmp_offset_read(cji) + is_key;
%                                 finalloc = (locval * fold) + (floor((tmp_offset_read(cji)-tmp_offset_read(1))/offset_inc)+1) + is_key;
%                                                                 else % regular sampled data ie angle gathers
%                                                                     finalloc = (locval * fold) + (floor((tmp_offset_read(cji)-tmp_offset_read(1))/job_meta.tkey_inc)+1) + is_key;
%                                                                 end
%                             else
                                finalloc = (locval * fold) + (floor((tmp_offset_read(cji)-tmp_offset_read(1))/offset_inc)+1) + is_key;
%                             end
                            if finalloc > (is_key+((n_traces_to_read)*fold)-1)
                                fprintf('mode of offset increment not correct, caused by very irregular fold or offset increment of %d\n',offset_inc);
                            elseif finalloc < 1
                                fprintf('mode of offset increment not correct, caused by very irregular fold or offset increment of %d\n',offset_inc);
                                finalloc = 1;
                                tmp_offset_read(cji) = cji;
                            end
                            traces{ii_block}(:,finalloc) = tmptracesblk(:,cji);
                            %ilxl_read{ii_block}(finalloc,:) = tmp_ilxl_read(cji,:);
                            offset_read{ii_block}(finalloc,:) = tmp_offset_read(cji);
                            lastoffval = tmp_offset_read(cji);
                            lastilxlval = sum(tmp_ilxl_read(cji,:));
                        end
                        
                        %replicate the list of inline and xlines to the
                        %fold of the gathers
                        for listiloop = 1:1:(listi-1) 
                            %needs to check the dimensions
                            ilxl_read{ii_block}(is_key+((listiloop-1)*fold):is_key+((listiloop)*fold)-1,:) = repmat(tmpilxllist(listiloop,:),fold,1);
                        end
                        
                        is_traces_to_read = (n_traces_to_read)*fold;
                    end
                    %alloffpres = size(offset_read{28},1)/fold
                    
                    s_key = s_key + n_traces_to_read;
                    i_counter = s_key; 
                    is_key = is_key + is_traces_to_read;
                    is_counter = is_key;
                                      
                end
            end
        end
    else
        % no traces
        seismic = segy_read_binary(strcat(job_meta.paths{1},blocks{ii_block}));
        %         if job_meta.is_gather == 1
        %             traces{ii_block} = zeros(seismic.n_samples,fold,'single');
        %             ilxl_read{ii_block} = zeros(fold,2,'int32');
        %             offset_read{ii_block} = int32(job_meta.tkey_min):int32(job_meta.tkey_inc):int32(job_meta.tkey_max)';
        %         else
        traces{ii_block} = zeros(seismic.n_samples,1,'single');
        ilxl_read{ii_block} = int32([0,0]);
        offset_read{ii_block} = int32(0);
        %         end
        
        % could to make zero entry
    end
end
%%
    if ztrunc == 1
        %truncate the cell of traces
        for ii_block = 1:1:size(blocks,1)
            traces{ii_block} = traces{ii_block}(zmin:zmax,:);
        end    
        % modify seismic structure to reflect the changes
        seismic.n_samples = (zmax - zmin )+ 1;
        seismic.firstsample = zmin; 
    end
%     if isfield(job_meta,'incl_poly')
%         for ii_block = 1:1:size(blocks,1)
%             il_r=ilxl_read(:,1);
%             xl_r=ilxl_read(:,2);
%             [IN ON] = inpolygon(coo(:,1),coo(:,2),il_poly,xl_poly);
%             INN=IN|ON;
%             
%             
%             
%             traces{ii_block} = traces{ii_block}(zmin:zmax,:);
%         end
%     end
    
    traces = cell2mat(traces);
    %cj fell over here
    ilxl_read = cell2mat(ilxl_read');
    zero_log = ilxl_read(:,1) ~= 0;
    ilxl_read = ilxl_read(zero_log,:);
    
    traces = traces(:,zero_log);
    
   % ilxl_read = ilxl_read';
    offset_read = cell2mat(offset_read');
    offset_read = offset_read(zero_log);
    offset_read = offset_read';
    
end
        
function [traces,ilxl_read,offset_read] = read_traces_segy(seismic,start_byte,n_traces_to_read)

    fid = fopen(char(seismic.filepath),'r','b');
    fseek(fid,start_byte,'bof');

    if seismic.file_type == 1
        % Convert traces from IBM32FP read as UINT32 into IEEE64FP (doubles) - need to make it singles
        %traces_tmp = fread(fid,[60+seismic.n_samples,n_traces_to_read],strcat(num2str(seismic.n_samples),'*uint32=>uint32'));
        traces_tmp = fread(fid,[60+seismic.n_samples,n_traces_to_read],'*uint32');
        
        
        %ilxl_read = traces_tmp(48:49,:)'; % what happens if the inline and crossline are not in this location  
        trchead = traces_tmp(1:60,:);
        [trace_header bytes_to_samples] = interpretbe(reshape(typecast(trchead(:),'uint16'),120,[]));
        
        % Inline / crossline compression step
        ilxl_read(:,1) = int32(trace_header(bytes_to_samples == seismic.pkey,:))';
        ilxl_read(:,2) = int32(trace_header(bytes_to_samples == seismic.skey,:))';
        offset_read = int32(trace_header(bytes_to_samples == seismic.tkey,:))';
        
        %offset_read = traces_tmp(10,:)';
        
        traces = single((1-2*double(bitget(traces_tmp(61:end,:),32))).*16.^ ...
        (double(bitshift(bitand(traces_tmp(61:end,:),uint32(hex2dec('7f000000'))),-24))-64).* ...
        (double(bitand(traces_tmp(61:end,:),uint32(hex2dec('00ffffff'))))/2^24));
    elseif seismic.file_type == 2 
        disp('This seismic file type is not currently supported. Please speak to Charles Jones.');
    elseif seismic.file_type == 5
        % Traces are IEEE32FP (singles)   
        traces = fread(fid,[60+seismic.n_samples,n_traces_to_read],strcat(num2str(seismic.n_samples),'*float32=>float32'));
        %trace_headers = typecast(single(reshape(traces(1:60,:),1,60*n_traces_to_read)),'int32');  
        %trace_headers = reshape(trace_headers,60,n_traces_to_read);
        
        trchead = typecast(single(reshape(traces(1:60,:),1,60*n_traces_to_read)),'uint16');  
        %trchead = reshape(trchead,120,n_traces_to_read);
        
        [trace_header bytes_to_samples] = interpretbe(reshape(trchead,120,[]));
        
        % Inline / crossline compression step
        ilxl_read(:,1) = int32(trace_header(bytes_to_samples == seismic.pkey,:))';
        ilxl_read(:,2) = int32(trace_header(bytes_to_samples == seismic.skey,:))';
        offset_read = int32(trace_header(bytes_to_samples == seismic.tkey,:))';
        
        %ilxl_read = trace_headers(48:49,:)';        
        %offset_read = trace_headers(10,:)';        
        traces = traces(61:end,:);            
    else
        disp('This seismic file type is not currently supported. Please speak to Charles Jones.');
    end

    fclose(fid);  
end


function [trace_header bytes_to_samples] = interpretbe(tmptrheader)
byte_type = [ ...
    2*ones(7,1); ones(4,1);
    2*ones(8,1); ones(2,1);
    2*ones(4,1); ones(46,1);
    2*ones(5,1); ones(2,1);
    2*ones(1,1); ones(5,1);
    2*ones(1,1); ones(1,1);
    2*ones(1,1); ones(2,1);
    2*ones(1,1); 2*ones(1,1)];

ntr = size(tmptrheader,2);
trace_header = zeros(91,ntr);
bytes_to_samples = zeros(91,1);

count =1;
for ii = 1:91
    bytes_to_samples(ii,1) = 2*count-1;
    if byte_type(ii) == 1
        trace_header(ii,:) = double(tmptrheader(count,:));
        count = count+1;
    elseif byte_type(ii) == 2
        trace_header(ii,:) = double(tmptrheader(count+1,:))*2^16 + double(tmptrheader(count,:)); % note this is big-endian and different to one in segy_make_structure
        count = count+2;
    end
end

trace_header(21,:) = trace_header(21,:)-2^16;

end