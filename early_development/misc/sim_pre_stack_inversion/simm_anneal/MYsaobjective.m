% Script to call the sort of objective function you wish to use
clear all
% Read data
% data = segy_read_files('/data/MAD/dtect/backup_brazil/TZA_seismics');
% trc1 = segy_read_traces(data{1},100,100,0,0);
% trc2 = segy_read_traces(data{2},100,100,0,0);
% trc3 = segy_read_traces(data{3},100,100,0,0);

% Estimate wavelets
w1 = estimate_zero_wavelet(trc1.data,data{1}.n_samples);
w2 = estimate_zero_wavelet(trc2.data,data{2}.n_samples);
w3 = estimate_zero_wavelet(trc3.data,data{3}.n_samples);

wavelets = horzcat(w1,w2,w3);
angles = [data{1}.angle*pi/180 data{2}.angle*pi/180 data{3}.angle*pi/180];
    
%N = 10; % discretisation in each direction
Sobs = vertcat(trc1.data,trc2.data,trc3.data);
avgIp = 4150;
avgIs = avgIp/2;
IpPRI = avgIp.*ones(length(trc1.data),1);
IsPRI = avgIs.*ones(length(trc1.data),1);
W1 = 0.8;
W2 = 0.1;
W3 = 0.1;

ObjectiveFunction = @(x) MYparameterized_objective(x,Sobs,IpPRI,IsPRI,W1,W2,W3,wavelets,angles);

% X0 is a column vector of length 2N
%X0 = [avgIp*ones(size(IpPRI)); avgIs*ones(size(IsPRI))];
X0 = [avgIp + 2000.*randn(size(IpPRI),1); avgIs + 1000.*randn(size(IsPRI),1)];
%r = avgIp + 2000.*randn(size(IpPRI),1);
% lb = avgIp
% ub =

options = saoptimset('Display','iter','MaxIter',5000);
[x,fval] = simulannealbnd(ObjectiveFunction,X0,[],[],options)
%[x,fval] = simulannealbnd(ObjectiveFunction,X0,lb,ub,options)
