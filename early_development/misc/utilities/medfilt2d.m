function [datafilt] = medfilt2d(data,xstepout,zstepout)

    [n_samples,n_traces] = size(data);
    idx = reshape((1:1:n_samples*n_traces),n_samples,n_traces);
    datafilt = zeros(n_samples,n_traces);
    
    % Zero pad with rows
    data = [zeros(zstepout,n_traces);data;zeros(zstepout,n_traces)];
    idx = [zeros(zstepout,n_traces);idx;zeros(zstepout,n_traces)];
    
    % Zero pad with columns
    data = [zeros(n_samples+2*zstepout,xstepout),data,zeros(n_samples+2*zstepout,xstepout)];
    idx = [zeros(n_samples+2*zstepout,xstepout),idx,zeros(n_samples+2*zstepout,xstepout)];
       
    for ii = 1:1+2*xstepout
        for jj = 1:1+2*zstepout
            
            n_row_blocks = floor((size(data,1)-(jj-1))/(1+2*zstepout));
            n_col_blocks = floor((size(data,2)-(ii-1))/(1+2*xstepout));
            
            datacell = mat2cell(data(jj:n_row_blocks*(1+2*zstepout)+(jj-1),ii:n_col_blocks*(1+2*xstepout)+(ii-1)),1+2*zstepout*ones(n_row_blocks,1),1+2*xstepout*ones(n_col_blocks,1));

            datacell = reshape(datacell,1,[]);

            datamat = sort(reshape(cell2mat(datacell),(1+2*xstepout)*(1+2*zstepout),[]));

            med = datamat(round(((1+2*xstepout)*(1+2*zstepout))/2),:);

            % Make logical mask of moving window size
            mask = zeros((1+2*zstepout),(1+2*xstepout));
            mask(round(((1+2*xstepout)*(1+2*zstepout))/2)) = 1;
            mask = repmat(mask,n_row_blocks+1,n_col_blocks+1);
            mask = [zeros(jj-1,size(mask,2));mask];
            mask = [zeros(size(mask,1),ii-1),mask];
            mask = mask(1:size(idx,1),1:size(idx,2));

            lin_idx = idx(logical(mask));

            datafilt(lin_idx(lin_idx~=0)) = med;
        end
        if ii == 1+2*xstepout
            fprintf('Done!\n');
        else
            fprintf('%d%%... ',round(100*(ii/(1+2*xstepout))));
        end
    end
end