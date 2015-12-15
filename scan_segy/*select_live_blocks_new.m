function select_live_blocks_new(job_meta_path)
% seismic_mat_path is the matlab .mat file
% find common traces in each volume
% report back byte locations of keys searched for

job_meta = load(job_meta_path);

i_block = 1;

if job_meta.is_gather == 0
    loop_index = size(job_meta.files,1);
else
    loop_index = size(job_meta.files,2);
end

i_counter = 1;

for i_block = 1:1:size(job_meta.files,1)
    seismic_mat_path = [job_meta.paths{1}, job_meta.files{i_block}];
    seismic = segy_read_binary(seismic_mat_path);

    for i_key = 1:1:size(job_meta.block_keys,1)
        expand_keys = expand_pst_keys(job_meta.block_keys(i_key,:),...
            job_meta.pkey_inc(1));
        
        is_live = lookup_block(seismic,expand_keys);
        
        if is_live > 0
            job_meta.liveblocks(i_counter) = i_key;
            i_counter = i_counter + 1;
        end
        
        
    end
end
 
% save(job_meta_path,'-struct','job_meta','-v7.3'); % Saves Seismic structure to mat file

end
      
function expand_keys = expand_pst_keys(keys,pkey_inc)
    %pkey_inc = mode(diff(seismic.trace_ilxl_bytes(:,1)));
    expand_pkeys = keys(1):pkey_inc:keys(2);

    %need to add key template for columns instead of numbering 5
    %skey_inc = mode(seismic.trace_ilxl_bytes(:,5));
    expand_skeys = keys(3):keys(4)-keys(3):keys(4);

    expand_pkeys = repmat(expand_pkeys,size(expand_skeys,2),1);
    expand_skeys = repmat(expand_skeys,size(expand_pkeys,2),1)';

    expand_keys = [expand_pkeys(:),expand_skeys(:)];
end

function is_live = lookup_block(seismic,expand_keys)
    is_live = 0;
    for i_key = 1:1:size(expand_keys,1)
        il_rows = seismic.trace_ilxl_bytes(:,1) == expand_keys(i_key,1);

        if sum(il_rows) > 0
            il_found = seismic.trace_ilxl_bytes(il_rows,:);
            xl_rows = (il_found(:,2) <= expand_keys(i_key,2)) & (il_found(:,4) >= expand_keys(i_key,2));

            if sum(xl_rows) > 0
                is_live = is_live + 1;
            else

                % no skey found

            end
        else

            % no peky found
        end
    end
end