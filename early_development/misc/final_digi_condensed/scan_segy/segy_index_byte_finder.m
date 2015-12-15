function vol_keys = segy_index_byte_finder(job_meta_path,keys,vol_index)
% seismic_mat_path is the matlab .mat file
% find common traces in each volume
% report back byte locations of keys searched for

job_meta = load(job_meta_path);
vol_index = str2double(vol_index);
% loop over all volumes
i_block = 1;
for i_files = 1:1:size(job_meta.files,1)
    if strfind(job_meta.files{i_files},job_meta.volumes{vol_index})
        seismic = segy_read_binary(strcat(job_meta.paths{1},job_meta.files{i_files}));

        expand_keys = expand_pst_keys(keys,...
            job_meta.pkey_inc(vol_index),job_meta.skey_inc(vol_index));

        keys_found = lookup_byte(seismic,expand_keys);
        if size(keys_found,2) == 3; % no entry if it does not occur in the block
            key_logic = keys_found(:,3) ~= 0;
            vol_keys{i_block,1} = keys_found(key_logic,:);
        else
            vol_keys{i_block,1} = 0;
        end
        i_block = i_block+1;
    end
end
    
% check if beyond range of survey
% check if not on grid
%end

end

function expand_keys = expand_pst_keys(keys,pkey_inc,skey_inc)
    %pkey_inc = mode(diff(seismic.trace_ilxl_bytes(:,1)));
    expand_pkeys = keys(1):pkey_inc:keys(2);

    %need to add key template for columns instead of numbering 5
    %skey_inc = mode(seismic.trace_ilxl_bytes(:,5));
    expand_skeys = keys(3):skey_inc:keys(4);

    expand_pkeys = repmat(expand_pkeys,size(expand_skeys,2),1);
    expand_skeys = repmat(expand_skeys,size(expand_pkeys,2),1)';

    expand_keys = [expand_pkeys(:),expand_skeys(:)];
end

function expand_keys = lookup_byte(seismic,expand_keys)
    %
    for i_key = 1:1:size(expand_keys,1)
        il_rows = seismic.trace_ilxl_bytes(:,1) == expand_keys(i_key,1);      

        if sum(il_rows) > 0
            il_found = seismic.trace_ilxl_bytes(il_rows,:);
            xl_rows = (il_found(:,2) <= expand_keys(i_key,2)) & (il_found(:,4) >= expand_keys(i_key,2));

            if sum(xl_rows) > 0 
                xl_found = il_found(xl_rows,:);

                % calculate number of traces away from starting byte
                [xl,ind] = min(xl_found(:,2));
                check_xl = (expand_keys(i_key,2) - xl)/xl_found(ind,5);
                
                bytes_per_sample = 4;
                trc_head = 240;
                trc_length = seismic.n_samples*bytes_per_sample;
                                
                n_traces_away = check_xl;
                n_bytes_away = n_traces_away*(trc_length+trc_head);
                
                expand_keys(i_key,3) = xl_found(ind,3)+n_bytes_away;                
            else
                % no skey found
                
            end
        else
            % no peky found
            












    end

    end
    

end
%     il_rows = seismic.trace_ilxl_bytes(:,1) == keys(:,1);    
% 
%     if sum(il_rows) > 0
%         il_found = seismic.trace_ilxl_bytes(il_rows,:);
% 
%    `xl_rows = (il_found(:,2) <= keys(:,2)) & (il_found(:,4) >= keys(:,2));
%     
%     if sum(xl_rows) > 0 
%         xl_found = il_found(xl_rows,:);
% 
%         offset_rows = (xl_found(:,6) <= keys(:,3)) & (xl_found(:,7) >= keys(:,3));
%         if sum(offset_rows) > 0
%             offset_found = xl_found(offset_rows,:);
%             
%             check_xl = (keys(:,2) - offset_found(:,2))/offset_found(:,5);
%             check_offset = (keys(:,3) - offset_found(:,6))/offset_found(:,8);
%             if check_xl/floor(check_xl) == 1 && check_offset/floor(check_offset) == 1
%                                 
%                 
%                 if seismic.file_type < 1 || seismic.file_type > 5 % need to break if not 1 or 5 because we don't handle it
%                     seismic.file_type = input('Non standard seismic file type. Please enter seismic file type (1 (IBM),2,3,4 or 5 (IEEE)): ', 's');
%                     seismic.file_type = str2num(seismic.file_type);
%                 else
%                     bytes_per_sample = 4;
%                 end
%                 trc_head = 240;
%                 trc_length = seismic.n_samples*bytes_per_sample;
%                 n_offsets = (offset_found(:,7)-offset_found(:,6))/offset_found(:,8);
%                 n_traces_away = (check_xl*(n_offsets+1))+check_offset;
%                 n_bytes_away = n_traces_away*(trc_length+trc_head);
%                 trace_byte = offset_found(:,3)+n_bytes_away;
%                 
%                 fid = fopen(char(seismic.filepath),'r','b');
%                 fseek(fid,trace_byte-240,'bof');
%                 
%                 traces_tmp = fread(fid,[60+seismic.n_samples,1],'*uint32');
%                 ilxl_read = traces_tmp(48:49,:)'; % what happens if the inline and crossline are not in this location        
%                 
%                 offset_read = traces_tmp(10,:)';
% 
%                 
%             else
%             
%             end
%             
%         else
%             
%         end     
%                 
%     else
%                    
%     end
% 
% end
%     
% end
% end