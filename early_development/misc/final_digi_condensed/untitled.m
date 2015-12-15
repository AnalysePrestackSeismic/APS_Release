        for i_row = 1:1:size(live_traces,2)            
            end_ind = end_ind + size(live_traces{i_row},2);
            %dataprob(live_traces{i_row},ii) = 1-datacdf(idxdataprob(start_ind:end_ind));      % lookup the probability scores from the cdf
            start_ind = end_ind + 1;
        end