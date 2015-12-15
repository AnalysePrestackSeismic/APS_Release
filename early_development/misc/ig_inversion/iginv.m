% Some setup parameters

% Scan seismic to get inline crossline byte locations
files_in = data_read('/segy/NOR/hegg_2010_geotrace/angle_stacks');

fprintf('\nThe following files are in directory: %s\n',files_in.path)
for i = 1:files_in.nfiles
    fprintf('File %d: %s\n',i,files_in.names{i});
end

index_files = input('Enter file numbers to process (in bracket []): ');
il_byte = 189;
xl_byte = 193;
aperture = 1; % full offset

% Scan each segy file (perhaps add flag to only scan 1 file)
for i = 1:length(index_files)
   seis{i} = segy_make_index(cell2mat(strcat(files_in.path,'/',files_in.names(index_files(i)))),il_byte,xl_byte);
end

ner_seis = segy_make_index('2318_angstk_raw_5_10_v2.sgy',189,193);
mid_seis = segy_make_index('2318_angstk_raw_10_15_v2.sgy',189,193);
far_seis = segy_make_index('2318_angstk_raw_15_20_v2.sgy',189,193);


% Make required positions for processing with an aperture of 1
ner_seis = segy_make_proc_pos(ner_seis,aperture);
mid_seis = segy_make_proc_pos(mid_seis,aperture);
far_seis = segy_make_proc_pos(far_seis,aperture);

% Invert in little groups based on proc
% Setup parameters for inversion that only happen once
nt = ner_seis.n_samples; % should eventually change this to not use nt at all 
taper = 10;
wrms = 1; 
tol = 1e-3;
iter = 500;
% Scales the RMS of the wavelet to be equal to wrms
% Setup angles
angles = [8,13,18];
nangle = length(angles);
n_pos = (2*aperture+1)^2;
angles_stepout = sort(repmat(angles,1,n_pos));
sin2theta = (sin(angles.*pi/180)).*(sin(angles.*pi/180));
sin2theta_stepout = (sin(angles_stepout.*pi/180)).*(sin(angles_stepout.*pi/180));

% Read group of traces
start_index = 1;

%vint = zeros(seismic.n_samples,seismic.n_traces);
%pos = zeros(2,seismic.n_traces);

%Iilf = zeros(ner_seis.n_samples,ner_seis.n_traces);
%Gilf = zeros(ner_seis.n_samples,ner_seis.n_traces);
inv_il = zeros(1,ner_seis.n_traces);
inv_xl = zeros(1,ner_seis.n_traces);
for ii=1:1:ner_seis.n_traces
   
    ner_traces = segy_read_traces(ner_seis,start_index,n_pos*ii,1);
    mid_traces = segy_read_traces(mid_seis,start_index,n_pos*ii,1);
    far_traces = segy_read_traces(far_seis,start_index,n_pos*ii,1); 

    start_index = start_index+n_pos; % this ensures that groups of n_pos are read
    
    fprintf('\nInverting trace %d (inline %d and crossline %d) of %d\n',...
        ii,...
        ner_traces.pos(1,(ner_traces.end_index-ner_traces.start_index)/2+1),...
        ner_traces.pos(2,(ner_traces.end_index-ner_traces.start_index)/2+1),...
        ner_seis.n_traces);

    % Apply tapering
    ner_traces.taper = ner_traces.data;
    mid_traces.taper = mid_traces.data;
    far_traces.taper = far_traces.data;

    for i=1:n_pos % need to check tapering is correct; looked weired when plotted?
        ner_traces.taper(end-taper:end,i) = (2/sqrt(pi))*ner_traces.data(end-taper:end,i).*cos(.5*pi*(0:taper)'/taper);
        ner_traces.taper(1:taper+1,i) = (2/sqrt(pi))*ner_traces.data(1:taper+1,i).*cos(.5*pi*(taper:-1:0)'/taper);
        mid_traces.taper(end-taper:end,i) = (2/sqrt(pi))*mid_traces.data(end-taper:end,i).*cos(.5*pi*(0:taper)'/taper);
        mid_traces.taper(1:taper+1,i) = (2/sqrt(pi))*mid_traces.data(1:taper+1,i).*cos(.5*pi*(taper:-1:0)'/taper);
        far_traces.taper(end-taper:end,i) = (2/sqrt(pi))*far_traces.data(end-taper:end,i).*cos(.5*pi*(0:taper)'/taper);
        far_traces.taper(1:taper+1,i) = (2/sqrt(pi))*far_traces.data(1:taper+1,i).*cos(.5*pi*(taper:-1:0)'/taper);
    end
    % End tapering

    S = [ner_traces.taper,mid_traces.taper,far_traces.taper];
    
    % Line fit to get start for Is
    Is = zeros(nt,1);
    Gs = Is;
    for i=1:nt
        tmp = robustfit(sin2theta_stepout,S(i,:),'ols');
        Is(i) = tmp(1,1);
        Gs(i) = tmp(2,1);
    end

    % Estimte wavelets (needs to be replaced)
    Sft = abs(fft(S));
    wraw = ifft(Sft/sqrt(nt),'symmetric');
    for i=1:n_pos*nangle
        wtmp = (2/sqrt(pi))*wraw(1:51,i).*cos(.5*pi*(0:50)'/(50.5));
        sn0 = wtmp(1)+2*sum(wtmp(2:end).*(-1).^(1:50)');
        wtmp(1) = wtmp(1)-sn0;
        wtmp = [wtmp(51:-1:2); wtmp];
        wtmp = wtmp/(norm(wtmp/wrms)/sqrt(length(wtmp)));
        w(:,i) = wtmp;
    end
    % End estimate wavelets

    % Make wavelet convolution matrix
    for i = 1:n_pos*nangle
        W_tmp(:,nt*(i-1)+1:nt*i) = convmtx(w(:,i),nt);
    end
    clip = (size(W_tmp,1)-nt)/2;
    W = W_tmp(clip+1:end-clip,:);
    W = sparse(W);
    % End make wavelet convolution matrix

    % Make whitenning matrix
    tmp_white = eye(nt);
    white = repmat(tmp_white,n_pos*nangle,2);
    lambdaI = 1;
    lambdaG = 1;
    white(:,1:nt) = white(:,1:nt)*lambdaI;
    white(:,1+nt:end) = white(:,1+nt:end)*lambdaG;
    white = sparse(white);
    % End make whitenning matrix

    % Make operator without linear trend
    op = [W',W'];
    for i=1:n_pos*nangle
        op(1+(i-1)*nt:i*nt,1+nt:end) = op(1+(i-1)*nt:i*nt,1+nt:end)*sin2theta_stepout(i);
    end
    op = op+white;
    % End Make operator without linear trend

    c = tan(16*pi/180);

    % Make operater with linear trend
    op_linfit = op;
    for i=1:n_pos*nangle
        op_linfit(1+(i-1)*nt:i*nt,1:nt) = op_linfit(1+(i-1)*nt:i*nt,1:nt)*(1-(sin2theta_stepout(i)/c));
    end

    % Tikhonov Regularization Matrix

    B = [-0.5,1,-0.5];
    B=repmat(B,2*nt,1);
    L = chol(spdiags(B,[-1,0,1],2*nt,2*nt));
    Iscale = 40;
    Gscale = 10;
    L(1:nt,:) = Iscale*L(1:nt,:);
    L(nt+1:end,:) = Gscale*L(nt+1:end,:);
    op_linfit = [op_linfit;L];

    [inv_linfit,~,residual,iter] = lsqr(op_linfit,[reshape(S,[],1);Is;zeros(nt,1)],tol,iter,[],[],[Is;zeros(nt,1)]);
    
    % need to save this out of the loop
    Iilf(:,ii) = inv_linfit(1:nt);
    Gilf(:,ii) = inv_linfit(nt+1:end)-(Iilf(:,ii)/c);
    
    % save position information
    inv_il(:,ii) = ner_traces.pos(1,(ner_traces.end_index-ner_traces.start_index)/2+1);
    inv_xl(:,ii) = ner_traces.pos(2,(ner_traces.end_index-ner_traces.start_index)/2+1);

end

% need to add write statement

fclose('all');




