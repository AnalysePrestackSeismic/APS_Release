
function [m,alpha_i] = conjugate_gradient(L,H,d,lambda,tol,N)

p = zeros(length(d),1);
m = zeros(length(d),1);
r = -d;

%[r1 c1] = size(L);
%[r2 c2] = size(d);
s_p = zeros(length(d),1);
s_m = zeros(length(d),1);
s_r = zeros(length(d),1);

%g_m = zeros(length(d),1);
%g_p = zeros(length(d),1);
%g_r = zeros(length(d),1);

    for n=1:1:N;
        g_m = L'*r - lambda*m;
        g_p = H'*g_m + lambda*p;
        g_m = H*g_p;
        g_r = L*g_m;
        rho = g_p'*g_p;

        if  n == 1;
            beta = 0;
            rho_0 = rho;
        else
            beta = rho/rho_hat;
        

            if beta < tol || rho/rho_0 < tol
                break
            end
        end
            s_p = g_p+beta*s_p;
            s_m = g_m+beta*s_m;
            s_r = g_r+beta*s_r;

            alpha = rho/(s_r'*s_r+lambda*((s_p'*s_p)-(s_m'*s_m)));
            

            p = p-alpha*s_p;
            m = m-alpha*s_m;
            r = r-alpha*s_r;    
            
            rho_hat = rho;
    end
end


