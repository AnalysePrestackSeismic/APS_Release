function [det_coef] = localcorr(data,zsmooth)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    [n_samples,n_traces] = size(data);
     
    S = (1/zsmooth)*spdiags(repmat([(1:1:zsmooth),(zsmooth-1:-1:1)],n_samples,1),[(-zsmooth+1:1:0),(1:1:zsmooth-1)],n_samples,n_samples);
    
    det_coef = zeros(size(data));
    
    for ii=1:n_traces-1
        
        a = data(:,ii);
        A = S*spdiags(a,0,n_samples,n_samples);
        b = data(:,ii+1);
        B = S*spdiags(b,0,n_samples,n_samples);

        det_coef(:,ii) =  ((A*b).*(B*a))./((A*a).*(B*b));     
    end

end

