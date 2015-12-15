function vrms_vint_dix_3d_(tvpairs)

z_step = 4; % 4 ms is 1 index step
min_z_out = 0;
max_z_out = 1240;
n_traces = size(tvpairs,2);
%for ii = 1:1:n_traces
ll_b = 0;
for ii = 1:1:n_traces
    ll = size(tvpairs{ii},1);
    if ll > ll_b        
        ll_b = ll;
        use_ii = ii;
    end
end

z_in = 1+(tvpairs{use_ii}(:,1)/z_step);
z_out = min_z_out:1:max_z_out;
n_out = length(z_out);

z_in_round = round(z_in);
z_in_round(end) = max_z_out;
n_in = length(z_in);
%end

for ii = 1:1:n_traces
    z_in = 1+(tvpairs{ii}(:,1)/z_step);
    tvpairs_int{ii}(:,2) = interp1(z_in,tvpairs{ii}(:,2),z_in_round,'spline');
    tvpairs_int{ii}(:,1) = z_in_round;   
end

tvpairs_int = cell2mat(tvpairs_int');
nsamp = n_in;
vint = vrms_vint_inv(tvpairs_int(:,1),tvpairs_int(:,2),n_traces,0.8,nsamp,n_traces)

end


function [vint] = vrms_vint_inv(z_in,vrms_in,n_traces_in,weight,nsamp,n_traces)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

data = z_in.*vrms_in.^2; % Data is t*Vrms^2 where k is sample point
%x0 = zeros(nsamp*n_traces,1);
data = [data;zeros(nsamp*n_traces,1);zeros(nsamp*n_traces,1)];
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
maxit = 200;
eta = 0.1;

%data = make_data(nsamp,nt,nxl,nil,double(nx),double(ny),double(nz),planarity,eta);

[vint,~,~,~] = cgls(@(x,flag)calculate_operator(x,flag,eta,z_in,nsamp,n_traces_in),data,0,tol,maxit,'true');    

imagesc(reshape(real(sqrt(vint(1:nsamp*n_traces))),nsamp,n_traces),[1.5 1.52])

end

function out = calculate_operator(x,flag,eta,z_in,nsamp,n_traces)
    %x = reshape(x,[],n_traces_in);    
    if flag == 1; % A*x 
        x = reshape(x(1:nsamp*n_traces),nsamp,n_traces);    
        % Perform integration
        ii_out = 2;
        i_trace = 1;
        %n_z_in = length(z_in);
        z_in(end+1) = 9999;
        while ii_out <= length(z_in);
            if z_in(ii_out-1) < z_in(ii_out)
                y(ii_out-1,1) = sum(x(1:z_in(ii_out-1),i_trace));
                ii_out = ii_out + 1;
            else
                i_trace = i_trace + 1;
                ii_out = ii_out + 1;
            end
        end
        %y = cumsum(x,1);
        x = x(:);
        % Perform temporal differentiation
        % Perform spatial differentiation
        dx = calculate_dx(nsamp,n_traces,1,'no_transp');
        dx = dx*x;
        dz = calculate_dz(nsamp,n_traces,1,'no_transp');
        dz = dz*x;
        %[dx,dz] = gradient(x,1);        
                   
        out = [y(:);eta.*dz(:);eta.*dx(:)];   
        
    elseif flag == 2; % A'*x
        x = reshape(x(1:nsamp*n_traces),nsamp,n_traces); 
        
        % To perform transpose
        x = flipud(x);
        ii_out = 2;
        i_trace = 1;
        %n_z_in = length(z_in);
        z_in(end+1) = 9999;
        while ii_out <= length(z_in);
            if z_in(ii_out-1) < z_in(ii_out)
                y(ii_out-1,1) = sum(x(1:z_in(ii_out-1),i_trace));
                ii_out = ii_out + 1;
            else
                i_trace = i_trace + 1;
                ii_out = ii_out + 1;
            end
        end        
        y = flipud(y);
        %y = flipud(cumsum(x,1));
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