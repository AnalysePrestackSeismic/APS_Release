function status = int_grad_inversion(angle_stacks_in,output_location,...
    proc_pos,tol,iter,taper,wrms)

    n_blocks = length(proc_pos);
    n_angles = length(angle_stacks_in);
    n_pos = (2*angle_stacks_in{1}.aperture+1)^2;
    
    angles_stepout = sort(repmat(angle_stacks_in.angle,1,n_pos));
    %sin2theta = (sin(angles.*pi/180)).*(sin(angles.*pi/180));
    sin2theta_stepout = (sin(angles_stepout.*pi/180)).*(sin(angles_stepout.*pi/180));
    
    start_index = 1;
    for i_block = 1:1:n_blocks
        % read in all traces in block
        traces{1:n_angles} = segy_read_traces(angle_stacks_in{1:n_angles},...
            proc_pos(i_block,4),proc_pos(i_block,5));
    
        % Apply tapering
        traces{1:n_files}.taper(end-taper:end) = ...
            (2/sqrt(pi))*traces{1:nfiles}.data(end-taper:end).*cos(.5*pi*(0:taper)/taper)';
        traces{1:n_files}.taper(1:taper+1) = ...
            (2/sqrt(pi))*traces{1:nfiles}.data(1:taper+1).*cos(.5*pi*(taper:-1:0)/taper)';
        
        n_groups = ntraces/n_pos;
        for ij=1:1:n_groups
           start = 1;
           for i=1:n_files
                S(:,start:n_pos*i) = traces{i}.taper;
                Sft(:,start:n_pos*i) = abs(fft(traces{i}.taper));
                start = start+n_pos;
           end
            
           start_index = start_index+n_pos;
        end
   
   start = 1;
   for i=1:n_files
        S(:,start:n_pos*i) = traces{i}.taper;
        Sft(:,start:n_pos*i) = abs(fft(traces{i}.taper));
        start = start+n_pos;
   end
   
   wraw = ifft(Sft/sqrt(nt),'symmetric');
   
    for i=1:n_pos*nangle
        wtmp = (2/sqrt(pi))*wraw(1:51,i).*cos(.5*pi*(0:50)'/(50.5));
        sn0 = wtmp(1)+2*sum(wtmp(2:end).*(-1).^(1:50)');
        wtmp(1) = wtmp(1)-sn0;
        wtmp = [wtmp(51:-1:2); wtmp];
        wtmp = wtmp/(norm(wtmp/wrms)/sqrt(length(wtmp)));
        if isnan(wtmp)
            w(:,i) = 0;
        else
            w(:,i) = wtmp;
        end
    end
    
     wangle = w*(spdiags(sin2theta_stepout',0,n_pos*nangle,n_pos*nangle));

 % Whiten wavelets
        lambdaI = 0.5;
        lambdaG = 0.5;
        ww = w;
        ww(51,:) = ww(51,:)+lambdaI;
        wanglew = wangle;
        wanglew(51,:) = wanglew(51,:)+lambdaG;

        % Add linear trend to wavelets
        c = tan(20*pi/180);
        wlinfitw = ww*(spdiags((1-(sin2theta_stepout'/c)),0,n_pos*nangle,n_pos*nangle));

        % Make wavelet convolution matri    %     Ww = sparse(Ww);
    %     Wanglew = sparse(Wanglew);
    %     Wlinfitw = sparse(Wlinfitw);
        % End make wavelet convolution matrix

        % Make operator without linear trend
    %     op = [Ww',Wanglew'];
        % End Make operator without linear trendces
        for i = 1:n_pos*nangle
    %         Ww_tmp(:,nt*(i-1)+1:nt*i) = convmtx(ww(:,i),nt);
            Wanglew_tmp(:,nt*(i-1)+1:nt*i) = convmtx(wanglew(:,i),nt);
            Wlinfitw_tmp(:,nt*(i-1)+1:nt*i) = convmtx(wlinfitw(:,i),nt);
        end
    %     Ww_tmp = sparse(Ww_tmp);
        Wanglew_tmp = sparse(Wanglew_tmp);
        Wlinfitw_tmp = sparse(Wlinfitw_tmp);
        clip = (size(Wanglew_tmp,1)-nt)/2;
    %     Ww = Ww_tmp(clip+1:end-clip,:);
        Wanglew = Wanglew_tmp(clip+1:end-clip,:);
        Wlinfitw = Wlinfitw_tmp(clip+1:end-clip,:);
    %     Ww = sparse(Ww);
    %     Wanglew = sparse(Wanglew);
    %     Wlinfitw = sparse(Wlinfitw);
        % End make wavelet convolution matrix

        % Make operator without linear trend
    %     op = [Ww',Wanglew'];
        % End Make operator without linear trend

        % Make operater with linear trend
        op_linfit = [Wlinfitw',Wanglew'];
        % End make operater with linear trend

        % Tikhonov Regularization Matrix
        B = [-0.5,1,-0.5];
        B = repmat(B,2*nt,1);
        L = chol(spdiags(B,[-1,0,1],2*nt,2*nt));
        Iscale = 10;
        Gscale = 10;
        L(1:nt,:) = Iscale*L(1:nt,:);
        L(nt+1:end,:) = Gscale*L(nt+1:end,:);
        op_linfit = [op_linfit;L];

        [inv_linfit,~,residual,iter] = lsqr(op_linfit,[reshape(S,[],1);zeros(nt,1);zeros(nt,1)],tol,iter,[],[],[traces{1}.data(:,(traces{1}.end_index-traces{1}.start_index)/2+1);zeros(nt,1)]);

        Iilf(:,i_block) = inv_linfit(1:nt);
        Gilf(:,i_block) = inv_linfit(nt+1:end)-(Iilf/c);
        
        IGilf.pos(1,ij) = traces.pos(1,(start_index+(ij*n_pos))/2); 
        IGilf.pos(2,ij) = traces.pos(2,(start_index+(ij*n_pos))/2);             
        start_index = start_index+n_pos;

   
    
    worker = num2str(proc_pos(i_block,2));
    block = num2str(proc_pos(i_block,3));
        
    % Write data
    data_I = strcat(output_location,'/','Iinv_worker_',worker,'_block_',block);
    file_I = fopen(data_I,'a');
    fwrite(file_I,Iilf,'float32');
    fclose(file_I);
    
    data_G = strcat(output_location,'/','Ginv_worker_',worker,'_block_',block);
    file_G = fopen(data_G,'a');
    fwrite(file_G,Gilf,'float32');
    fclose(file_I);

    % Write position information
    meta_out = strcat(output_location,'/','vrms_pos_worker_',worker,'_block_',block);
    dlmwrite(meta_out,vint.pos','delimiter','\t');
        
    end

end