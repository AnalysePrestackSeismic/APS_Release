function [ mask ] = low_amp_mute( data_tmpb )
%LOW_AMP_MUTE Summary of this function goes here
%   Detailed explanation goes here

%old route
%data_tmp = data_tmpb(:);
% Find zones where data is zero (due to mute angle mute functions)
% data_zeros = abs(data_tmp) < abs(mean(data_tmp)/range(data_tmp)); 
%
%     data_tmp = abs(vol_traces(:,:,ii));
%     data_tmp = data_tmp(:);
%     mu3 = mean(data_tmp); % Data mean
%     sigma3 = std(data_tmp); % Data standard deviation
%     outliers = (data_tmp - mu3) > 2*sigma3;
%
% data_tmpb = vol_traces(:,:,ii);
% ==================================================================
% this applies a short vertical smoother to the data and then compares
% those values to the mean of each constant time, the mean is moothed over
% several time steps to get a bigger average, this give a %variation to the
% running mean

filtw = [1 2 2 2 1]/8;
data_tmpb(isnan(data_tmpb)) = 0;
data_tmpb = abs(data_tmpb);
%data_tmpb(data_tmpb == 0) = NaN;
%cjkl = nanmean(abs(cjk),2); can get a better mean with nanmean, but low is
%good for this mute anyway
data_tmpa = conv2(filtw,1,data_tmpb,'same');

tmpmean = mean(data_tmpb,2);
cjk = repmat(tmpmean,1,size(data_tmpb,2));
data_tmpb(bsxfun(@ge,data_tmpb,(tmpmean*1.5)) == 1) = cjk(bsxfun(@ge,data_tmpb,(tmpmean*1.5)) == 1);
finmean = mean(data_tmpb,2);

filt_smo =  ones(1,21)/21;
tmpmean_smo = conv(finmean,filt_smo,'same');
finmean(finmean > tmpmean_smo) = tmpmean_smo(finmean > tmpmean_smo);
tmpd = bsxfun(@rdivide,data_tmpa,finmean);
%
% ==================================================================
% now apply a short vertival smoother and pick the point of low difference
%filtww = [1 2 3 2 1]/9;
tmpf = conv2(filtw,1,abs(tmpd),'same') >  0.05;

%tmpe = conv2(filtww,1,abs(tmpd),'same');
%data_tmpb(tmpe < 0.1) = 0;
%[FX,FY] = gradient(tmpd);
%imagesc(abs(FY))
% percentage change from the pre filtered version
%tmpe = abs((data_tmpa  ./ (abs(data_tmpb) - data_tmpa))  + 1);

% ==================================================================
% 
% low pick the left hand side and the right hand side mutes
[~,lhs] = max(single(tmpf'));
lhs(end+1:1:end+length(filt_smo)) = lhs(end);
lhs = conv(lhs,filt_smo,'same');
lhs = floor(lhs(1:size(data_tmpb,1)))-1;
lhs(lhs < 0) = 0;

[~,rhs] = max(single(tmpf(:,end:-1:1)'));
rhs = size(data_tmpb,2)-rhs;
rhs(end+1:1:end+length(filt_smo)) = rhs(end);
rhs = conv(rhs,filt_smo,'same');
rhs = floor(rhs(1:size(data_tmpb,1))+0.5)+1;
rhs(rhs > size(tmpf,2)) = size(tmpf,2);

% 2d filtering
%s = [1 2 1; 0 0 0; -1 -2 -1];
%s = ones(5,5)/25;
%data_tmpa = conv2(data_tmpb,s,'same');

%outline = union(lhs, rhs, 'rows');
%mask = 1:size(data_tmpb,1)*size(data_tmpb,2);
%mask = reshape(mask,size(data_tmpb,2),[])';

% ==================================================================
% now use the mutes to make a mask
% use expansion in a for loop
% mask = zeros(size(data_tmpb,1),size(data_tmpb,2),'single');
% for ij = 1: size(data_tmpb,1)
%     mask(ij,lhs(ij):rhs(ij)) = 1;
% end
% 
% make the mask using bsxfun rather than a for loop
mask = ones(size(data_tmpb,1),size(data_tmpb,2),'single');
mask = cumsum(mask,2);
mask(bsxfun(@ge,mask,lhs') == 0)  = 0;
mask(bsxfun(@le,mask,rhs') == 0)  = 0;
%mask(mask ~= 0) = 1;
%
% add in points for a blank data to the mask, ie the hard zeros round the
% edge
mask = data_tmpa .* mask;

%make mask 0's and 1's
mask(mask ~= 0) = 1;

% % apply the mutes to blank low amp data
% data_tmpb(reshape(data_zeros,[],size(data_tmpb,2))) = 0;
% data_tmpb = data_tmpb .* mask;

end

