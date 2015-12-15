function [wb_idx] = water_bottom_flatten_lite(traces)
    exponent = 50;
    inverse_gain = (size(traces,1):-1:1)'.^exponent;

    traces = bsxfun(@times,traces,inverse_gain);

    [~, wb_idx] = max(traces);
end

