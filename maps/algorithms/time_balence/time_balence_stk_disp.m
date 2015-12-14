function trim_data_filt = time_balence_stk( trim_data )
%% Defination: TIME_BALENCE scales a gather or setion from the envelope of the amplitudes
% Input:  
% Output:
% Writes to Disk:

%to make the envelope tend to the value +/- 2000
%   Detailed explanation goes here
%%
filt_smo =  ones(1,3)/3;
filt_smo_big =  ones(1,91)/91;
%filttraces = [1 2 2 3 3 3 2 2 1]/19;
trim_data_filt = trim_data;
scalto = 126;
scaltoneg = -126;
toplim = 127;
toplimneg = -127;

for ckk = 1:size(trim_data,2)
    
% find the max of the data across the gather and smooth
td_max = trim_data(:,ckk);
%td_max =  conv(td_max,filttraces,'same');

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

%use the peaks logical to get the values and indexes to then interpolate to
%make and envelope which includes the signal, but preserves the wiggles in
%the dataset
itpsval = td_max(diffsign < 0);
itpslocs = single(1:size(trim_data,1))';
itpslocsin = itpslocs(diffsign < 0);

% interpolate to make the envelope only using fast linear interp
posenv = double(interp1q(itpslocsin,itpsval,itpslocs));


%use the mins logical to get the values and indexes to then interpolate to
%make and envelope which includes the signal, but preserves the wiggles in
%the dataset
mitpsval = td_max(mdiffsign > 0);
mitpslocs = single(1:size(trim_data,1))';
mitpslocsin = mitpslocs(mdiffsign > 0);

% interpolate to make the min envelope only using fast linear interp
mposenv = double(interp1q(mitpslocsin,mitpsval,mitpslocs));

maxenv = max([abs(posenv) abs(mposenv)],[],2);
maxenv = conv(maxenv,filt_smo_big,'same');
maxenv((size(maxenv,1)-size(filt_smo_big,2)+1):end,1) = maxenv((size(maxenv,1)-size(filt_smo_big,2)+1),1);
% now make the scaler to make envlope all fit the value 2000
scalepos = scalto ./ maxenv;

%scalepos = bsxfun(@rdivide,2000,posenv);
scalepos(isnan(scalepos)) = 0;

% apply a median filter to remove sudden jumps in scaling and the small
% averaging to make sure it is a smooth scalar
%scalepos = medfilt3nt(scalepos,15,0);
%scalepos =  conv(scalepos,filt_smo_big,'same');


%Apply the scaling to the input data
%trim_data_filt(ckk) = bsxfun(@times,scalepos,trim_data(ckk));
trim_data_filt(:,ckk) = trim_data(:,ckk).*scalepos;
trim_data_filt((trim_data_filt(:,ckk) > scalto),ckk) = toplim;
trim_data_filt((trim_data_filt(:,ckk) < scaltoneg),ckk) = toplimneg;
end

end

