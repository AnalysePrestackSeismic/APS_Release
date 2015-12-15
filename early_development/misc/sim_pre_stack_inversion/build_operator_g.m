function [G] = build_operator_g(wavelets,diff,c1,c2,c3)
    %UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here

    % wavelet convolution matrix
    cropw = ceil((lw-1)/2); % 
    nrw = nt+lw-1; %
        for i=1:1:nw
            tmp = convmtx(wavelets(:,i),nt);
            wavelets_conv{i} = sparse(tmp(cropw:nrw-cropw,:)); 
            g_tmp{i} = [(c1(i,1)*wavelets_conv{i}*diff) (c2(i,1)*wavelets_conv{i}*diff) (c3(i,1)*wavelets_conv{i}*diff)];       
        end

    % operator or gradient matrix
    G = vertcat(g_tmp{1:end}); % gradient matrix for inversion problem

end