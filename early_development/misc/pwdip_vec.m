function [dip] = pwdip_vec(data)
% Uses plane wave prediction filters to solve for dip. This method solves
% the nonlinear inverse problem cast by Sergy Fomel in 'Applications of
% plane-wave destruction filters' to find the dip which optimally predicts
% the next trace from the previous. In this formulation the trace is
% modelled as a plane wave, allowing dip to be handled by addition to the
% phase term.

    [n_samples,n_traces] = size(data);
    
    dip = zeros(n_samples,n_traces-1);
    
    % Nuisance parameters. Lambda controls smoothness of the dip estimate.
    % Prewhite stabilizes noisy areas. iter_max is the maximum iteration of
    % the nonlinear inversion.
    lambda = 10000;
    prewhite = 100;
    iter_max = 10;
           
    d = data(:);

    for ii = 1:iter_max

        [C, Cd] = pwf(dip(:,jj),n_samples);

        % Update the model paramters (dip) by find the least squares
        % best fit for dm=(inv((Cd*D)+prewhite))*(-C*d) under the
        % constraint that lambda*diff(dm)=0.
        dm = (spdiags([(Cd*d)+prewhite*ones(n_samples,1),lambda*ones(n_samples,1),-lambda*ones(n_samples,1)],[0,n_samples-1,n_samples],n_samples,2*n_samples)'\[-C*d;zeros(n_samples,1)]);

        dip(:,jj) = dip(:,jj)+dm;
    end
    dip = 2*dip;
end
        
function [C, Cd] = pwf(dip,n_samples)
% Calculate filter coefficients for plane wave destruction and the
% associated filter of derivatives needed for the nonlinear inversion.

    c11 = -((1+dip).*(2+dip))/12;
    c12 = ((1-dip).*(2-dip))/12;
    c21 = -((2+dip).*(2-dip))/6;
    c22 = ((2+dip).*(2-dip))/6;
    c31 = -((1-dip).*(2-dip))/12;
    c32 = ((1+dip).*(2+dip))/12;
    
    % Plane wave destruction filter
    C = [spdiags([c11(:),c21(:),c31(:)],[-1,0,1],n_traces*(n_samples-1),n_traces*(n_samples-1)),spdiags([c12(:),c22(:),c32(:)],[-1,0,1],n_traces*(n_samples-1),n_traces*(n_samples-1))];
    
    cd11 = ((-2*dip)-3)/12;
    cd12 = ((2*dip)-3)/12;
    cd21 = dip/3;
    cd22 = -dip/3;
    cd31 = (3-(2*dip))/12;
    cd32 = ((2*dip)+3)/12;
    
    % Filter of derivatives
    Cd = [spdiags([cd11(:),cd21(:),cd31(:)],[-1,0,1],n_traces*(n_samples-1),n_traces*(n_samples-1)),spdiags([cd12(:),cd22(:),cd32(:)],[-1,0,1],n_traces*(n_samples-1),n_traces*(n_samples-1))];

end