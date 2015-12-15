clear all
%clc

file_in = input('Enter filename: ');
nil = input('Enter number of inlines: ');
nxl = input('Enter number of xlines: ');
nt = input('Enter number of time samples: ');
nwin = input('Enter number of time samples for window length: ');
olap = input('Enter time window overlap (0 = none, 0.1 = 10%, etc): ');

fprintf('\nReading data...\n');

ns = nil*nxl*nt;
ntrace = ns/nt;
fid = fopen(file_in);
data = fread(fid,[nt,ntrace],'float32');
fclose(fid);

%file_in = 'vol_grad_raw_wheeler_flat.bin';
%file_in = 'vol_int_raw_seabed_flat.bin';

fprintf('\nCorrecting undefineds...\n');

data = data';
for k=1:ntrace
    if nansum(abs(data(k,:))) == 0
        data(k,:) = NaN;
    else
        for i=1:nt
            if data(k,i)>1e29
                data(k,i) = NaN;
            end
        end
    end
    if single((ntrace/50)*floor(k/(ntrace/50))) == k;
        fprintf('%d / %d = %.0f%%\n',k,ntrace,(k*100)/(ntrace));
    end
end

win_shift = nwin-(floor(nwin*olap));
nloop = ceil((nt/win_shift)-1);
add_col = ((nloop+1)*win_shift)-nt;
average = zeros(1,nloop);
variance = zeros(1,nloop);

data = [data,NaN*zeros(ntrace,add_col)];
count = 1;

fprintf('\nCalculating statistics...\n');

for k = 0:win_shift:((nloop-1)*win_shift)
    data_tmp = reshape(data(:,k+1:k+nwin),[],1);
    average(1,count) = nanmean(data_tmp);
    variance(1,count) = nanvar(data_tmp); 
    stddev = sqrt(variance);
    count = count+1;
end

clearvars data_tmp

fprintf('\nInterpolating...\n');

x = nwin/2:win_shift:nloop*win_shift;
xi = 1:1:nt;
averagei = (interp1(x,average,xi,'pchip','extrap'));
variancei = (interp1(x,variance,xi,'pchip','extrap'))';
stddevi = (interp1(x,stddev,xi,'pchip','extrap'));

% fprintf('\nSaving variance to file...\n');

%variancei = variancei'*ones(1,ntrace);
% fid1 = fopen(sprintf('variance_win%d_olap%d_%s',nwin,nwin-win_shift,file_in),'a');
% for k = 1:ntrace
%     fwrite(fid1,variancei,'float32');
% end
% fclose(fid1);

data = data(:,1:end-add_col);

% fprintf('\nCalculating Q-Q matrix...\n');
% 
% p = (1:1:99);
% qq = zeros(nt,length(p));
% normd = randn(1000000,1);
% 
% for k = 1:nt
%     qq(k,:) = prctile(data(:,k),p)-prctile(averagei(1,k)+stddevi(1,k)*normd,p);
% end

fprintf('\nClearing memory...\n');

clearvars -except averagei stddevi data nwin win_shift file_in ntrace

fprintf('\nCalculating cdf...\n');

prob = zeros(size(data));
for k = 1:ntrace
    prob(k,:) = (2*cdf('Normal',abs(data(k,:)-averagei)+averagei,averagei,stddevi))-1;
    if single((ntrace/50)*floor(k/(ntrace/50))) == k;
            fprintf('%d / %d = %.0f%%\n',k,ntrace,(k*100)/(ntrace));
    end
end

clearvars data

fprintf('\nSaving average and cdf to files...\n');

%averagei = averagei'*ones(1,ntrace);
prob = prob';

% fid2 = fopen(sprintf('average_win%d_olap%d_%s',nwin,nwin-win_shift,file_in),'a');
% for k=1:ntrace
%     fwrite(fid2,averagei,'float32');
% end
% fclose(fid2);

fid3 = fopen(sprintf('cdf_win%d_olap%d_%s',nwin,nwin-win_shift,file_in),'w');
fwrite(fid3,prob,'float32');
fclose(fid3);

fprintf('\nComplete\n');