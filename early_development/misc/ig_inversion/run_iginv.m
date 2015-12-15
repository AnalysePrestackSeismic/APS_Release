clear all

seis = segy_read_files('/lustrecache/TZA/dtect/kussini_2012_full_seq_pgs/Seismics/angle_stacks');

aperture = 0; 
n_files = length(seis);
hor_flag = input('Do you want to use a horizon? [1 - Yes, 0 - No]: ');
for i = 1:n_files
    if hor_flag == 1 && i > 1
        seis{i}.horpath = seis{1}.horpath;
        seis{i}.time_start = seis{1}.time_start;
        seis{i}.proc_samples = seis{1}.proc_samples;
    end
    seis{i} = segy_make_proc_pos(seis{i},aperture,hor_flag);   
    angles(i) = seis{i}.angle;
estart_index = start_index+n_pos;nd

% Start parameters for inversion
nangle = length(angles); % number of angle stacks to process
n_pos = (2*aperture+1)^2; % number of positions to process per trace
if hor_flag
    nt = seis{1}.proc_samples;
else
    nt = seis{1}.n_samples;
end
taper = 10; % tapering length for cosine taper at start and end of trace
wrms = 1; % scaling term
tol = 1e-3; % tolerance of residual
iter = 500; % max. number of iterations

angles_stepout = sort(repmat(angles,1,n_pos));
sin2theta = (sin(angles.*pi/180)).*(sin(angles.*pi/180));
sin2theta_stepout = (sin(angles_stepout.*pi/180)).*(sin(angles_stepout.*pi/180));

% Error check on input trace numbers needed

% Read group of traces
start_index = 1;

for ii=1:1:seis{1}.n_traces
    
    zero_check = 0;
    for ij = 1:n_files
        traces{ij} = segy_read_traces(seis{ij},start_index,n_pos*ii,1,hor_flag);
        if traces{ij}.data == 0
            zero_check = zero_check + 1;
        end
    end
    
    start_index = start_index+n_pos; % this ensures that groups of n_pos are read
    
    if zero_check ~= n_files
    
         % Apply tapering to data
        for i=1:n_files
            % traces{1:n_files}.taper = traces{1:n_files}.data; need to check
            % how to logical index a cell array
            traces{i}.taper = traces{i}.data;
        end

        for i=1:n_files
            traces{i}.taper(end-taper:end) = (2/sqrt(pi))*traces{i}.data(end-taper:end).*cos(.5*pi*(0:taper)/taper)';
            traces{i}.taper(1:taper+1) = (2/sqrt(pi))*traces{i}.data(1:taper+1).*cos(.5*pi*(taper:-1:0)/taper)';
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

        % Make wavelet convolution matrices
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

        Iilf = inv_linfit(1:nt);
        Gilf = inv_linfit(nt+1:end)-(Iilf/c);
        
    else
        
        Iilf = zeros(nt,1);
        Gilf = zeros(nt,1);
        
    end

        file_i = fopen('inv_intercept.bin','a');

        fwrite(file_i,Iilf,'float32');

        file_g = fopen('inv_gradient.bin','a');
        fwrite(file_g,Gilf,'float32');

        fclose all; 

        if (floor(ii/100) == ii/100)
            fprintf('Processing trace %d of %d (%d%% complete)\n',ii,seis{1}.n_traces,round((ii/seis{1}.n_traces)*100));
        end

    %     inv_il(:,ii) = ner_traces.pos(1,(ner_traces.end_index-ner_traces.start_index)/2+1);
    %     inv_xl(:,ii) = ner_traces.pos(2,(ner_traces.end_index-ner_traces.start_index)/2+1);
    
end
    
fclose all; 




