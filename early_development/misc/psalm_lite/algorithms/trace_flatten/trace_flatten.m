function [traces_flat] = trace_flatten(traces,wb_idx)

    for kk = 1:length(wb_idx)
        traces_flat(:,kk) = circshift(traces{ii}(:,kk),-wb_idx(kk));
        traces_flat(end-wb_idx(kk):end,kk) = 0;
    end

end