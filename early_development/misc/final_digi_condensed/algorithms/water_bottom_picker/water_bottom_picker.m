function [wb_idx] = water_bottom_picker(traces,add_to_scan)
% Water bottom flatten function
% Input:
% Output:   

    exponent = 50;
    inverse_gain = (size(traces,1):-1:1)'.^exponent;
    traces = bsxfun(@times,traces,inverse_gain);
    [~, wb_idx] = max(traces);
    wb_idx = wb_idx-5;
    %wb_idx = medfilt1(wb_idx,20);
    % add header to scan file
    if add_to_scan == 1
        
    else
        
    end
end

