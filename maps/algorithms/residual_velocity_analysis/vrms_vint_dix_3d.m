function vrms_vint_dix_3d(tvpairs)

z_step = 4; % 4 ms is 1 index step
min_z_out = 0;
max_z_out = 1240;
n_traces = size(tvpairs,2);
i_col = 1;
for ii = 1:1:n_traces    
    z_in = 1+(tvpairs{ii}(:,1)/z_step); % convert to sample index
    z_in_round = floor(z_in); % round down to nearest sample
    z_in_round(end) = max_z_out;
    tvpairs{ii}(:,2) = interp1(z_in,tvpairs{ii}(:,2),z_in_round,'spline'); % interpolate vels onto new axis  
    tvpairs{ii}(:,1) = z_in_round;    
    tvpairs{ii}(:,3) = i_col;
    i_col = i_col + 1;
end
z_out = min_z_out:1:max_z_out;
nsamp = length(z_out);

tvpairs = cell2mat(tvpairs');

% restriction matrix indexes
m = zeros(nsamp,n_traces); % model matrix



vint = vrms_vint_inv(tvpairs(:,1),tvpairs(:,2),n_traces,0.8,nsamp,n_traces)

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
R = sub2ind([nsamp,n_traces],tvpairs(:,1)',tvpairs(:,3)'); % make linear indices to select
%data = make_data(nsamp,nt,nxl,nil,double(nx),double(ny),double(nz),planarity,eta);

[vint,~,~,~] = cgls(@(x,flag)calculate_operator(x,flag,z_in,n_traces_in),data,0,tol,maxit,'true');    

end

function out = calculate_operator(x,flag,z_in,n_traces_in)
    %x = reshape(x,[],n_traces_in);    
    if flag == 1; % A*x 
        x = reshape(x,[],n_traces_in);    
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
    
        % Perform temporal differentiation
        % Perform spatial differentiation
        [dx,dz] = gradient(x,1);        
                   
        out = [y(R);dz(:);dx(:)];
        
    elseif flag == 2; % A'*x
        % Perform integration
        z_in = flipud(z_in); % reverse z_in
        ii_out = 2;
        i_trace = 1;
        %n_z_in = length(z_in);
        z_in(end+1) = 9999;
        while ii_out <= length(z_in);
            if z_in(ii_out-1) > z_in(ii_out)
                y(ii_out-1,1) = sum(x(1:z_in(ii_out-1),i_trace));
                ii_out = ii_out + 1;
            else
                i_trace = i_trace + 1;
                ii_out = ii_out + 1;
            end
        end
        
    end
end