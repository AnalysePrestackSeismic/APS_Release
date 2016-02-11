function [vint] = vrms_vint_dix_3d_grid(gridtv,eta)
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.
%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
%   Detailed explanation goes here

[nsamp,n_traces] = size(gridtv);

z_in = 1:1:nsamp;
%data = gridtv.^2;

% for i_trace = 1:1:n_traces
%     data(:,i_trace) = z_in'.*data(:,i_trace);
% end
data = bsxfun(@times,gridtv.^2,z_in'); % Data is t*Vrms^2 where k is sample point
x0 = zeros(nsamp*n_traces,1);
data = [data(:);zeros(nsamp*n_traces,1);zeros(nsamp*n_traces,1)];
%x0 = 
% for iz = 1:n_in
%     C(iz,1:z_in_round(iz)) = ones(1,z_in_round(iz));
% end

%D = spdiags([weight*ones(max_z_out,1) -2*weight*ones(max_z_out,1) weight*ones(max_z_out,1)],[-1 0 1],max_z_out,max_z_out);
%D = spdiags([weight*-1*ones(max_z_out,1) weight*ones(max_z_out,1)],[0 1],max_z_out,max_z_out);

%vint = real(sqrt(lsqr([C;D],[data;zeros(max_z_out,1)])));

%plot(sqrt(m),-1*(0:1:max_z_out-1))

disp('Calculate vint')
tol = 1e-4;
maxit = 10;
%eta = 0.5;

%data = make_data(nsamp,nt,nxl,nil,double(nx),double(ny),double(nz),planarity,eta);

[vint,~,~,~] = cgls(@(x,flag)calculate_operator(x,flag,eta,nsamp,n_traces),data,0,tol,maxit,'true');    
imagesc(reshape(real(sqrt(vint(1:nsamp*n_traces))),nsamp,n_traces),[1.5 1.52])
 
end

function out = calculate_operator(x,flag,eta,nsamp,n_traces)
    %x = reshape(x,[],n_traces_in);    
    if flag == 1; % A*x 
        x = reshape(x(1:nsamp*n_traces),nsamp,n_traces);    
        % Perform integration
        y = cumsum(x,1);
        x = x(:);
        % Perform temporal differentiation
        % Perform spatial differentiation
        dx = calculate_dx(nsamp,n_traces,1,'no_transp');
        dx = dx*x;
        dz = calculate_dz(nsamp,n_traces,1,'no_transp');
        dz = dz*x;
        %[dx,dz] = gradient(x,1);        
        % perform restriction           
        out = [y(:);eta.*dz(:);eta.*dx(:)];   
        
    elseif flag == 2; % A'*x
        x = reshape(x(1:nsamp*n_traces),nsamp,n_traces); 
        
        % To perform transpose
        x = flipud(x);
        y = flipud(cumsum(x,1));
        x = x(:);
        % Perform temporal differentiation
        % Perform spatial differentiation
        dx = calculate_dx(nsamp,n_traces,1,'transp');
        dx = dx*x;
        dz = calculate_dz(nsamp,n_traces,1,'transp');
        dz = dz*x;
        %[dx,dz] = gradient(x,1);        
        out = [y(:);eta.*dz(:);eta.*dx(:)];           
      
 
    end
end

function dx = calculate_dx(nt,nxl,nil,transp_flag) % x crossline
    %data should be a column vector 
    %from 3d volume with dimensions nt,nxl,nil
   if strcmp(transp_flag,'no_transp')
        % dx
        nsamp = nt*nxl*nil; 
        dx1 = [-1*ones(nt,1); zeros((nxl-2)*nt,1); ones(nt,1)];
        dx2 = [ones(nt,1); 1/2*ones((nxl-2)*nt,1); zeros(nt,1)];
        dx2 = [ones(nt,1); dx2(1:end-nt,1)];        
        dx3 = [-1/2*ones((nxl-2)*nt,1); -1*ones(nt*2,1)];          
       
        dx = spdiags([dx3,dx1,dx2],[-nt 0 nt],nt*nxl,nt*nxl);
        
        for i_loop = 1:nil;
            dx_cell{i_loop} = dx;
        end
                
        dx = blkdiag(dx_cell{:});
        %dx = dx*data;
        
        %s_dx= dx*double(vol_3d_v);
        %s_dx = reshape(s_dx,nt,nxl,nil);        
   elseif strcmp(transp_flag,'transp')
        % dx' 
        nsamp = nt*nxl*nil;
        dx1 = [-1*ones(nt,1); zeros((nxl-2)*nt,1); ones(nt,1)];
        dx2 = [ones(nt,1); 1/2*ones((nxl-2)*nt,1); zeros(nt,1)];
        %dx2 = [ones(nt,1); dx2(1:end-nt,1)];        
        dx3 = [-1/2*ones((nxl-2)*nt,1); -1*ones(nt*2,1)];    
        dx = spdiags([dx2,dx1,dx3],[-nt 0 nt],nt*nxl,nt*nxl);
        
        for i_loop = 1:nil;
            dx_cell{i_loop} = dx;
        end
                
        dx = blkdiag(dx_cell{:});
        %dx = dx*data;
   end
end

function dy = calculate_dy(nt,nxl,nil,transp_flag) % y inline
    % data should be a column vector 
    % from 3d volume with dimensions nt,nxl,nil
    if strcmp(transp_flag,'no_transp')
        nsamp = nt*nxl*nil; 
        dy1 = [-1*ones(nt*nxl,1); zeros(nsamp-2*nt*nxl,1); ones(nt*nxl,1)];
        dy2 = [ones(2*nt*nxl,1); 1/2*ones(nt*nxl*nil-2*nt*nxl,1);];
        dy3 = [-1/2*ones(nt*nxl*nil-2*nt*nxl,1); -1*ones(2*nt*nxl,1)];       
      
        dy = spdiags([dy3,dy1,dy2],[-nxl*nt 0 nxl*nt],nt*nxl*nil,nt*nxl*nil);             
        %dy = dy*data;
    elseif strcmp(transp_flag,'transp')
        nsamp = nt*nxl*nil; 
        dy1 = [-1*ones(nt*nxl,1); zeros(nsamp-2*nt*nxl,1); ones(nt*nxl,1)];
        dy2 = [ones(nt*nxl,1); 1/2*ones(nt*nxl*nil-nt*nxl,1);];
        dy3 = [-1/2*ones(nt*nxl*nil-nt*nxl,1); -1*ones(nt*nxl,1)];         
      
        dy = spdiags([dy2,dy1,dy3],[-nxl*nt 0 nxl*nt],nt*nxl*nil,nt*nxl*nil);             
        %dy = dy*data;    
    end
end

function dz = calculate_dz(nt,nxl,nil,transp_flag)
        % data should be a column vector 
        % from 3d volume with dimensions nt,nxl,nil
        if strcmp(transp_flag,'no_transp')
            % dz
            nsamp = nt*nxl*nil;       
            dz1 = repmat([-1 zeros(1,nt-2) 1],1,nsamp/nt)';
            dz2 = repmat([-1/2*ones(1,nt-2) -1 0],1,nsamp/nt)';
            dz3 = repmat([1 1/2*ones(1,nt-2) 0],1,nsamp/nt)';
            dz3 = [1;dz3(1:end-1)];
            dz = spdiags([dz2,dz1,dz3],[-1 0 1],nsamp,nsamp);
            %dz = dz*data;
        elseif strcmp(transp_flag,'transp')
            % dz'
            nsamp = nt*nxl*nil;
            dz1 = repmat([-1 zeros(1,nt-2) 1],1,nsamp/nt)';
            dz2 = repmat([-1/2*ones(1,nt-2) -1 0],1,nsamp/nt)';
            dz3 = repmat([1 1/2*ones(1,nt-2) 0],1,nsamp/nt)';       
            dz2 =  [-1/2;dz2(1:end-1)];
            dz = spdiags([dz3,dz1,dz2],[-1 0 1],nsamp,nsamp);
            %dz = dz*data;
        end      
end