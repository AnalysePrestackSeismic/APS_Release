% Version 0 - MCMC Trace by Trace

% OUT  w0:        estimated wavelet (normalised with standard deviation of noise)
%      vn0:       estimated average noise variance
%      w0_phase:  estimated phase in degrees
%      data_out:  deconvolved data
%      qc_kurtosis:  kurtosis values for all angles; size(n_angle,2)
%                      format: (angle, kurtosis)

clear
close all
clc

% read data
file_in = input('''Input file name'' (segy) = ');
endian = input('Endian format (''b'' = big, ''l'' = little) = ');
[Data,SegyTraceHeaders,SegyHeader]=ReadSegyConstantTraceLength(file_in,'endian',endian);
% file_out_rhat = input('''Wiener Output file name'' (always big endian on output) = ');
file_out_rm = input('''MCMC Output file name'' (always big endian on output) = ');

% nt = number of time samples
% nrec = number of traces
[nt, nrec] = size(Data);
dt = SegyHeader.dt/1.e6;

% parameters needed by bliss_kpe in order to find wavelet from data
tw = input('Wavelet time length (s) = ');
nw_h_in = round(tw/dt/2);
get_phase_in = input('get phase of wavelet (1 = yes, 0 = no) = ');
damp_in = input('damping constant in Wiener filtering (>1) = ');
% bliss_kpe function call
[w0,vn0,w0_phase,data_out,qc_kurtosis]=bliss_kpe(Data,nw_h_in,get_phase_in,damp_in);

% Input:
% - x: the data, assumed to obey the model x = conv(w,r) + noise
% - w: the wavelet (filter), assumed known and of short length
% - sn: the noise standard deviation, assumed known
% - n1: number of warm-up cycles in MCMC
% - n2: number of cycles in MCMC to compute the deconvolution result
% - verbose: set this variable is to true to get some printout
% Output:
% - rhat: Wiener deconvolution in the time domain (a column vector with
%	  length = length(x) + length(w)-1)
% - rm: MCMC deconvolution, same size as rhat
% * if n1, hence n2, are not given, only rhat is computed.

% mcmc decon function call
% need to loop over each trace and calculate the deconvolved trace
% w0 from bliss_kpe
% 
warm_in = input('number of warm-up cycles = ');
cycle_in = input('number of cycles in MCMC to compute decon. = ');

w0 = w0*sqrt(vn0);

for irec=1:nrec
%for irec=10:20
    [rhat(:,irec),rm(:,irec)] = mcmcdecon(Data(:,irec),w0,sqrt(vn0),warm_in,cycle_in);  
end

% plot original data
figure(1)
imagesc(1:nrec,(1:nt)*dt,Data); colormap(flipud(gray))
title('Original Data')
xlabel('CDP')
ylabel('Time (s)')

% plot Wiener data
figure(2)
imagesc(1:nrec,(1:nt)*dt,rhat); colormap(flipud(gray))
title('Wiener deconvolution')
xlabel('CDP')
ylabel('Time (s)')

% plot MCMC decon. data
figure(3)
imagesc(1:nrec,(1:nt)*dt,rm); colormap(flipud(gray))
title('MCMC non-linear deconvolution')
xlabel('CDP')
ylabel('Time (s)')

% save results
% note Segy on output always big endian
% WriteSegyStructure(file_out_rhat,SegyHeader,SegyTraceHeaders,rhat);
WriteSegyStructure(file_out_rm,SegyHeader,SegyTraceHeaders,rm);

rm_shift = circshift(rm,-nw_h_in);
rm_shift = rm_shift(1:nt,1:nrec);

save results.mat rm rm_shift w0 w0_phase warm_in cycle_in

fid = fopen('mcmc_decon.txt','wt');
fprintf(fid,'Wiener output file name: %s\n',file_out_rhat);
fprintf(fid,'MCMC Output file name: %s\n',file_out_rm);
fprintf(fid,'Warm-up cycles: %d\n',warm_in);
fprintf(fid,'MCMC cycles: %d\n',cycle_in);
fprintf(fid,'estimated average noise variance: %.8g\n',vn0);
fprintf(fid,'half length of wavelet: %d\n',nw_h_in);
fprintf(fid,'get phase of wavelet (1 = yes, 0 = no): %d\n',get_phase_in);
fprintf(fid,'damping constant in Wiener filtering: %d\n',damp_in);
fprintf(fid,'Wavelet: \n');
fprintf(fid,'%.8g\n',w0);
fclose(fid);