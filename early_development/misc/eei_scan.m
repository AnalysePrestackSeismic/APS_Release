function [sce] = eei_scan(data,nil,nxl,nt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[nfiles,~] = size(data);
if nfiles ~= 2
    error('Data input must be a cell array of size 2x1, with cell (1,1) being an intercept volume and cell (2,1) being a gradient volume')
end

I = reshape(data{1,1},[],1);
G = reshape(data{2,1},[],1);

I = (I-mean(I))/std(I);
G = (G-mean(G))/std(G);

[~,indexi] = sort(I.^2);
[~,indexg] = sort(G.^2);

for i = 1:length(I)
    if (indexi(i,1) >= round(0.01*length(I))) || (indexg(i,1) >= round(0.01*length(G)))
        Ibg(i,1) = NaN;
        Gbg(i,1) = NaN;
    else
        Ibg(i,1) = I(i,1);
        Gbg(i,1) = G(i,1);
    end
end

chi = (-90:5:90)*(pi/180);
stdev = zeros(length(chi),1);

for k = 1:length(chi)
    eei = Ibg*cos(chi(k))+Gbg*sin(chi(k));
    stdev(k,1) = nanstd(eei);
end

% for k = 1:length(chi)
%     eei = I*cos(chi(k))+G*sin(chi(k));
%     eei = sort(abs(eei));
%     diff_mean(k,1) = mean(eei(round(0.75*length(eei))+1:end))/mean(eei(1:round(0.75*length(eei))))));
% end

[~,index] = max(stdev);
sce = chi(index)*180/pi;

end
