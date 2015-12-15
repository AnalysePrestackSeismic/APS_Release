function [m] = cg_sreg(L,H,d,lambda,tol,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


p = zeros(size(L,1),1);
m = zeros(size(L,2),1);
r = -d;

for ii = 1:N
    gm = L'*r - lambda*m;
    gp = H'*gm + lambda*p;
    gm = H*gp;
    gr = L*gm;
    rho = gp'*gp;
    
    if ii == 1
        beta = 0;
        rho_hat = rho;
    else
        beta = rho./rho_hat;
    end
    
    if or(beta<tol,rho./rho_hat<tol)
        break
    end
    
    sp = gp + beta*sp;
    sm = gm + beta*sm;
    sr = gr + beta*sr;
    
    alpha = rho/(sr'*sr + lambda*(sp'*sp - sm'*sm));
    
    p = p - alpha*sp;
    m = m - alpha*sm;
    r = r - alpha*sr;
    
    rho_hat = rho;
end

