function [wb_idx] = water_bottom_picker(traces)
% Water bottom flatten function
% Input:
% Output:   


exponent = 50;
    inverse_gain = (size(traces,1):-1:1)'.^exponent;

    traces = bsxfun(@times,traces,inverse_gain);

    [~, wb_idx] = max(traces);
end

