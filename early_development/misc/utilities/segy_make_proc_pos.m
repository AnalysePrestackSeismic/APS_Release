function seismic = segy_make_proc_pos(seismic,aperture)
% for very large files this may need to run in a distributed manner    

    nfiles = length(seismic);
     
    for ii=1:1:nfiles
        
        if isfield(seismic{ii},'process')
            seismic{ii} = rmfield(seismic{ii},'process');
        end            
        
        if aperture == 0
            seismic{ii}.process = seismic{ii}.trace_pointers;
            seismic{ii}.aperture = aperture;
        else
            seismic{ii}.aperture = aperture;
            n_pos = (2*aperture+1)^2;

            pos_diff_xl = repmat(aperture:-1:-aperture,sqrt(n_pos),1);
            pos_diff_il = pos_diff_xl';
            pos_diff_il =  pos_diff_il(:).*seismic{ii}.il_inc; % turn into vector
            pos_diff_xl =  pos_diff_xl(:).*seismic{ii}.xl_inc;

            % do without for loop
            for ij=1:1:seismic{ii}.n_traces
                seismic{ii}.process_il(:,ij) = seismic{ii}.trace_pointers(ij,1)+pos_diff_il;
                seismic{ii}.process_xl(:,ij) = seismic{ii}.trace_pointers(ij,2)+pos_diff_xl;  
            end 

            seismic{ii}.process_il = seismic{ii}.process_il(:);
            seismic{ii}.process_xl = seismic{ii}.process_xl(:);
            seismic{ii}.process(:,1) = seismic{ii}.process_il;
            seismic{ii}.process(:,2) = seismic{ii}.process_xl;

            seismic{ii} = rmfield(seismic{ii}, 'process_il');
            seismic{ii} = rmfield(seismic{ii}, 'process_xl');

            [Lia, Locb] = ...
            ismember(seismic{ii}.process(:,1:2),seismic{ii}.trace_pointers(:,1:2), 'rows');

            seismic{ii}.process(:,3) = NaN(length(seismic{ii}.process(:,1)),1);
            seismic{ii}.process(Lia,3) = seismic{ii}.trace_pointers(nonzeros(Locb),3);

        end
        
        seismic{ii}.npos = length(seismic{ii}.process(:,1));

    end   

end