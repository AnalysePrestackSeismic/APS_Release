function [dip] = pwd_2D_dip_example(data,zsmooth,xsmooth)
% Estimates dip from a section using inversion of plane wave destruction
% (pwd) filters, L: L(-dip)*t1-L(dip)*t2 = 0. The linearized problem is: 
% Gm = d, where G = diag(dL(-dip)*t1-dL(dip)*t2) and 
% d = -(L(-dip)*t1-L(dip)*t2). 'dip' in the linearized problem is the dip
% estimate for the current iteration and m is the update needed to add to
% dip for the next iteration. dL is the derivative of L with respect to
% dip.

% INPUTS
% data = seismic section
% zsmooth = regularization weight controlling vertical smoothness
% xsmooth = regularization weight controlling horizontal smoothness

% OUTPUTS
% dip = seismic dip (units depend upon with data is in time or depth)

    [n_samples,n_il] = size(data);   
    dip = zeros(n_samples,n_il);   
    itermax = 4;
    
    % Iterativly solve the linearized inverse problem for incremental dip
    % update
    for kk = 1:itermax
        d = zeros(n_samples,n_il);
        G = zeros(n_samples,n_il);
        
        % Rebuild pwd filters from current dip estimate and update d and G
        for ii = 1:n_il-1
            [L, dL] = pwf(dip(:,ii));
            left = [data(1,ii);data(:,ii);data(end,ii)];
            right = [data(1,ii+1);data(:,ii+1);data(end,ii+1)];
            d(:,ii) = (L(:,3).*left(1:end-2,:)+L(:,2).*left(2:end-1,:)+L(:,1).*left(3:end,:))-...
                (L(:,1).*right(1:end-2,:)+L(:,2).*right(2:end-1,:)+L(:,3).*right(3:end,:));
            G(:,ii) = (dL(:,3).*left(1:end-2,:)+dL(:,2).*left(2:end-1,:)+dL(:,1).*left(3:end,:))-...
                (dL(:,1).*right(1:end-2,:)+dL(:,2).*right(2:end-1,:)+dL(:,3).*right(3:end,:));
        end
        d = -d(:);
        G = G(:);
        if kk == 1
            % Calculate energy of r.h.s. for normalizing regularization
            % weights
            energy = norm(G.*d);
        end
        % Solve inversion for m and add it to current dip estimate
        [m,~] = pcg(@(x)operate(x,G,zsmooth*energy,xsmooth*energy,n_samples,n_il),G.*d,1e-1,100);
        dip = dip+reshape(m,n_samples,n_il);
    end
    dip = 2*dip;
end

function y = operate(x,G,lambdaz,lambdax,n_samples,n_il)
% Performs [diag(G)'*diag(G)+lambda*D'*D]*m, where D = 2D differentation
% matrix

    % [diag(G)'*diag(G)]*m
    y = G.*G.*x; % n.b. diag(G)*x == G.*x and diag(G)' == diag(G)
    
    % [lambda*D'*D]*m
    x = reshape(x,n_samples,n_il);
    zedges(1:2,:) = cat(1,x(1,:),x(end,:));
    xedges(:,1:2) = cat(2,x(:,1),x(:,end));
    y = y+lambdaz*reshape(-diff(cat(1,x,zedges(2,:)),1,1)+diff(cat(1,zedges(1,:),x),1,1),[],1);
    y = y+lambdax*reshape(-diff(cat(2,x,xedges(:,2)),1,2)+diff(cat(2,xedges(:,1),x),1,2),[],1);
end

function [L,dL] = pwf(dip)
% Calculate pwd filter coefficients and the associated filter of
% derivatives 
    
    % Plane wave destruction filter
    L = [((1-dip).*(2-dip))/12,((2+dip).*(2-dip))/6,((1+dip).*(2+dip))/12];
    % Filter of derivatives
    dL = [(2*dip-3)/12,-dip/3,(2*dip+3)/12];
end