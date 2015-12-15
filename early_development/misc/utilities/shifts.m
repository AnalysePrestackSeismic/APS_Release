function [dip] = shifts(seis)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[seis_dewhorl] = dewhorl(seis);

real_seis = real(hilbert(seis_dewhorl));
imag_seis = imag(hilbert(seis_dewhorl));




end


function [seis_dewhorl] = dewhorl(seis)

    % Get the dimensions of the seismic input
    [ns ntr] = size(seis);
    
    whorl_idx = zeros(ns,1);
    
    zero_blocks = 0.5*(1+cumsum(abs(diff([0;sign(seis(:,ii))]))));
    extrema_blocks = 0.5*(1+cumsum(abs(diff([0;sign(diff([0;seis(:,ii)]))]))));
    
    for jj = 1:max(zero_blocks)
        zero_block_idx = ismember(zero_blocks,jj);
        extrema_blocks_in_zero_block = extrema_blocks(zero_block_idx);
        if extrema_blocks_in_zero_block(end)-extrema_blocks_in_zero_block(end)>2
            whorl_idx(zero_block_idx) = 1;
        end
    end
    
    
    
%     for ii=1:ntr
%     
%         % Create matrix that increments each time a zero is crossed
%         zero_blocks = cumsum(abs(diff([0;sign(seis(:,ii))])));
% 
%         % Create matrix that increments each time an extrema is crossed
%         extrema_blocks = cumsum(abs(diff([0;sign(diff([0;seis(:,ii)]))])));
% 
%         % Create matrix that locates the edges of extrema blocks
%         extrema_edges = logical(abs(conv2(sign(diff([0;seis(:,ii)])),[-1;1],'same')));
% 
%         % Create matrix that locates the starting edges of zero crossing blocks
%         start_zero_edges = logical(abs(conv2(sign(seis(:,ii)),[0;-1;1],'same')));
% 
%         % Create matrix that locates the ending edges of zero crossing blocks
%         end_zero_edges = logical(abs(conv2(sign(seis(:,ii)),[-1;1],'same')));
% 
%         % Create matrix that locates the zero crossing blocks that contain an whorl (a whorl exists when extrema blocks increment more than once within a single zero crossing block)
%         whorl_zero_blocks = zero_blocks(start_zero_edges);
%         end_extrema_blocks = extrema_blocks(end_zero_edges);
%         start_extrema_blocks = extrema_blocks(start_zero_edges);
%         if end_extrema_blocks(1) < start_extrema_blocks(1)
%             end_extrema_blocks = end_extrema_blocks(2:end);
%         end
%         if length(start_extrema_blocks) > length(end_extrema_blocks)
%             start_extrema_blocks = start_extrema_blocks(1:end-1);
%         end
%         whorl_zero_blocks = whorl_zero_blocks(end_extrema_blocks-start_extrema_blocks>2);
% 
%         % Find the indicies of the zero crossing blocks that contain a whorl
%         whorl_zero_idx = ismember(zero_blocks,whorl_zero_blocks);
%         
%         whorl_zero_edges = logical(diff([0;whorl_zero_idx]));
% 
%         % Find all the extremas within zero crossing blocks that contain a whorl
%         whorl_extrema_edges = logical(extrema_edges.*whorl_zero_idx);
% 
%         % Make a linear index of sample number for later interpolation
%         z_idx = (1:1:ns)';
% 
%         % Find the z indicies of all the extremas within zero crossing blocks that contain a whorl
%         whorl_extrema_edges_z_idx = z_idx(whorl_extrema_edges);
%         whorl_zero_edges_z_idx = z_idx(whorl_zero_edges);
%         [~,nearest_whorl_extrema_edge] = min(abs(bsxfun(@minus,whorl_extrema_edges_z_idx,whorl_zero_edges_z_idx')));
%         whorl_extrema_edges_z_idx_keep = whorl_extrema_edges_z_idx(unique(nearest_whorl_extrema_edge));
% 
% %         % Make a mask to select only the first and last z indicies of extremas within zero crossing blocks that contain a whorl
% %         whorl_extrema_edges_z_idx_keep_mask = logical(repmat([1;0;1],ceil(size(whorl_extrema_edges_z_idx,1)/3),1));
% %         whorl_extrema_edges_z_idx_keep_mask = whorl_extrema_edges_z_idx_keep_mask(1:size(whorl_extrema_edges_z_idx,1),:);
% % 
% %         % Select only the first and last z indicies of extremas within zero crossing blocks that contain a whorl
% %         whorl_extrema_edges_z_idx_keep = whorl_extrema_edges_z_idx(whorl_extrema_edges_z_idx_keep_mask);
% 
%         % Make a logical index table of all the samples to keep, excluding the whorls
%         keep_z_idx = or(ismember(z_idx,whorl_extrema_edges_z_idx_keep),...
%             floor(0.5*cumsum(ismember(z_idx,whorl_extrema_edges_z_idx_keep))) == 0.5*cumsum(ismember(z_idx,whorl_extrema_edges_z_idx_keep)));
% 
%         % Select only the z indicies of all points not in a whorl
%         z_idx_dec = z_idx(~keep_z_idx);
% 
%         % Select only the seismic samples of all points not in a whorl
%         seis_dec = seis(~keep_z_idx,ii);
% 
%         % Interpolate seismic samples across whorls
%         seis_dewhorl(:,ii) = interp1(z_idx_dec,seis_dec,z_idx,'pchip',0);
%     end

end
    
 
