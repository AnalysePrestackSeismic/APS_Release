function trim_data_filt = time_balence_stk( trim_data )
%% Defination: TIME_BALENCE scales a gather or setion from the envelope of the amplitudes
% Input:  
% Output:
% Writes to Disk:
%scalepeak = single(32700);
%scalepeak = single(23000);
scalepeak = single(110);
%to make the envelope tend to the value +/- 2000
%   Detailed explanation goes here
%%
filt_smo =  ones(1,3)/3;
%filttraces = [1 2 2 3 3 3 2 2 1]/19;
%trim_data_filt = trim_data;
trim_data_filt = zeros(size(trim_data,1),size(trim_data,2),'int8');
itpslocs = single(1:size(trim_data,1))';
itpslocstest = itpslocs;
origcount = itpslocs;

for ckk = 1:size(trim_data,2)
    
% reset some variables
itpslocs = origcount;
itpslocstest = origcount;
% find the max of the data across the gather and smooth
td_max = trim_data(:,ckk);
%td_max =  conv(td_max,filttraces,'same');

%figure(2); plot(td_max);
%hold all

% find the first and second derivatives of the max
max_1st_deriv = diff(td_max);
max_2nd_deriv = diff(td_max,2);

% apply a signum filter to get samples at zero crossings and make only 1
% and 0's
sign_1_deriv = sign(max_1st_deriv);
sign_1_deriv(sign_1_deriv == 0) = 1;

% find the point where sign of 1st deriv changes
diffsign = diff(sign_1_deriv);
mdiffsign = diffsign;

% set the point to zero where second derivative is positive and pad, this
% finds the peaks in the max dataset
diffsign(sign(max_2nd_deriv) > 0) = 0;
diffsign = [1;diffsign];

% set the point to zero where second derivative is positive and pad, this
% finds the mins in the max dataset
mdiffsign(sign(max_2nd_deriv) <= 0) = 0;
mdiffsign = [1;mdiffsign];
orig_mdiffsign = mdiffsign;
%use the peaks logical to get the values and indexes to then interpolate to
%make and envelope which includes the signal, but preserves the wiggles in
%the dataset
itpsval = td_max(diffsign < 0);
itpslocsin = itpslocs(diffsign < 0);
%itpsval = single(scalepeak ./ itpsval);

% interpolate to make the envelope only using fast linear interp
posenv = double(interp1q(itpslocsin,itpsval,origcount));
%plot(posenv)

%==========================================================================
% calc the envelope of the envelope
max_1st_deriv = diff(posenv);
max_2nd_deriv = diff(posenv,2);
sign_1_deriv = sign(max_1st_deriv);
sign_1_deriv(sign_1_deriv == 0) = 1;
diffsign = diff(sign_1_deriv);
mdiffsign = diffsign;
diffsign(sign(max_2nd_deriv) > 0) = 0;
diffsign = [1;diffsign];
mdiffsign(sign(max_2nd_deriv) <= 0) = 0;
mdiffsign = [1;mdiffsign];
itpsval = posenv(diffsign < 0);
itpslocsin = itpslocs(diffsign < 0);
% add small white noise to avoid division errors/mess
itpsval = itpsval + (trimmean(itpsval,10)*0.001);
% interpolate to make the envelope only using fast linear interp
posenv = double(interp1q(itpslocsin,itpsval,origcount));
%plot(posenv)
%==========================================================================

%==========================================================================
% % calc the envelope of the envelope of the envelope
% max_1st_deriv = diff(posenv);
% max_2nd_deriv = diff(posenv,2);
% sign_1_deriv = sign(max_1st_deriv);
% sign_1_deriv(sign_1_deriv == 0) = 1;
% diffsign = diff(sign_1_deriv);
% mdiffsign = diffsign;
% diffsign(sign(max_2nd_deriv) > 0) = 0;
% diffsign = [1;diffsign];
% mdiffsign(sign(max_2nd_deriv) <= 0) = 0;
% mdiffsign = [1;mdiffsign];
% itpsval = posenv(diffsign < 0);
% itpslocsin = itpslocs(diffsign < 0);
% % add small white noise to avoid division errors/mess
% itpsval = itpsval + (trimmean(itpsval,10)*0.001);
% % interpolate to make the envelope only using fast linear interp
% posenv = double(interp1q(itpslocsin,itpsval,origcount));
%plot(posenv)
%==========================================================================

itpsval = single(scalepeak ./ itpsval);
%==========================================================================
%use the mins logical to get the values and indexes to then interpolate to
%make and envelope which includes the signal, but preserves the wiggles in
%the dataset
mitpsval = td_max(orig_mdiffsign > 0);
%mitpslocs = single(1:size(trim_data,1))';
mitpslocsin = itpslocs(orig_mdiffsign > 0);

% interpolate to make the min envelope only using fast linear interp
mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
%plot(mposenv)

%==========================================================================
% calc the envelope of the envelope
max_1st_deriv = diff(mposenv);
max_2nd_deriv = diff(mposenv,2);
sign_1_deriv = sign(max_1st_deriv);
sign_1_deriv(sign_1_deriv == 0) = 1;
diffsign = diff(sign_1_deriv);
mdiffsign = diffsign;
diffsign(sign(max_2nd_deriv) > 0) = 0;
diffsign = [1;diffsign];
mdiffsign(sign(max_2nd_deriv) <= 0) = 0;
mdiffsign = [1;mdiffsign];
mitpsval = td_max(mdiffsign > 0);
mitpslocsin = itpslocs(mdiffsign > 0);
% add small white noise to avoid division errors/mess
mitpsval = mitpsval + (trimmean(mitpsval,10)*0.001);
% interpolate to make the envelope only using fast linear interp
mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
%plot(mposenv)
%==========================================================================

% % calc the envelope of the envelope of the envelope
% max_1st_deriv = diff(mposenv);
% max_2nd_deriv = diff(mposenv,2);
% sign_1_deriv = sign(max_1st_deriv);
% sign_1_deriv(sign_1_deriv == 0) = 1;
% diffsign = diff(sign_1_deriv);
% mdiffsign = diffsign;
% diffsign(sign(max_2nd_deriv) > 0) = 0;
% diffsign = [1;diffsign];
% mdiffsign(sign(max_2nd_deriv) <= 0) = 0;
% mdiffsign = [1;mdiffsign];
% mitpsval = td_max(mdiffsign > 0);
% mitpslocsin = itpslocs(mdiffsign > 0);
% % add small white noise to avoid division errors/mess
% mitpsval = mitpsval + (trimmean(mitpsval,10)*0.001);
% % interpolate to make the envelope only using fast linear interp
% mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
%plot(mposenv)
%==========================================================================

mitpsval = single(-scalepeak ./ mitpsval);

%==========================================================================

itpslocs(mitpslocsin) = mitpsval;
itpslocs(itpslocsin) = itpsval;
alllocsin = itpslocstest(itpslocs ~= itpslocstest);
allvalsin = itpslocs(itpslocs ~= itpslocstest);

allvalsin(isnan(allvalsin)) = 0;


% interpolate to make the min envelope only using fast linear interp
scalepos = single(interp1q(alllocsin,allvalsin,itpslocstest));
%figure(3); plot(scalepos)
%maxenv = max([abs(posenv) abs(mposenv)],[],2);


% now make the scaler to make envlope all fit the value 2000
%scalepos = 2000 ./ maxenv;
%scalepos = bsxfun(@rdivide,2000,posenv);
scalepos(isnan(scalepos)) = 0;
scalepos_restore = single(1./scalepos);

scalemean = trimmean(scalepos,50);

% apply a median filter to remove sudden jumps in scaling and the small
% averaging to make sure it is a smooth scalar
%scalepos = medfilt3nt(scalepos,15,0);
%scalepos =  conv(scalepos,filt_smo,'same');


%Apply the scaling to the input data
%trim_data_filt(ckk) = bsxfun(@times,scalepos,trim_data(ckk));
trim_data_filt(:,ckk) = int8(trim_data(:,ckk).*scalepos);
%figure(4); plot(trim_data_filt(:,ckk));

diffout = single(single(trim_data_filt(:,ckk)).*scalepos_restore);
uncerror =  (trim_data(:,ckk) - diffout) ./ trim_data(:,ckk);
uncerror(isnan(uncerror)) = 0;
uncerror(scalepos > scalemean*20) = 0;
unerrcount = sum(uncerror ~= 0);
avgpercterr = (sum(abs(uncerror)))/unerrcount;
%fprintf('%-10.8f percent error in reconstruction of dataset\n',avgpercterr);
fprintf('%-10.8f percent error in reconstruction of dataset, one in %d, using %d scalars\n',avgpercterr,round(100/avgpercterr),size(alllocsin,1));

%figure(5); plot(diffout);
%hold all
%plot(trim_data(:,ckk));

end
%     td_max = max(trim_data,[],2);
%     td_max =  conv(td_max,filt_smo,'same');
% 
%     cjdiff2 = diff(td_max,2);
%     cjdiff = diff(td_max);
%     cjdiffb = sign(cjdiff);
%     cjdiffb(cjdiffb == 0) = 1;
%     diffsign = diff(sign(cjdiff));
%     cjdiff2sign = sign(cjdiff2);
%     diffsign(cjdiff2sign > 0) = 0;
%     diffsign = [1;diffsign];
%     %td_max(diffsign == 0) = 0;
%     itpsval = td_max(diffsign < 0);
%     itpslocs = single(1:size(trim_data,1))';
%     itpslocsin = itpslocs(diffsign < 0);
%     posenv = double(interp1q(itpslocsin,itpsval,itpslocs));
%     %scalepos = 2000 ./ posenv;
%     scalepos = bsxfun(@rdivide,2000,posenv);
%     scalepos(isnan(scalepos)) = 0;
%     scalepos =  conv(scalepos,filttraces,'same');
%     scaltd = bsxfun(@times,scalepos,trim_data);

end

