function outzone = stk_event_zones( trim_data , percent_to_drop, passzonewidth)
%TIME_BALENCE scales a gather or setion from the envelope of the amplitudes
%to make the envelope tend to the value +/- 2000
%   Detailed explanation goes here

filt_smo =  ones(1,3)/3;
passfilt = ones(1,passzonewidth)/passzonewidth;
%filttraces = [1 2 3 3 3 3 3 2 1]/21;
trim_data_filt = trim_data;

% a gaussian filter
t = linspace(-1,1,25)';
filttraces = zeros(length(t),1);
a = 1;
filttraces = sqrt(pi)/a*exp(-(pi*t/a).^2);
filttraces = filttraces/sum(filttraces);


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
    
    % try to select zones of events
    residmean = min(maxenv(maxenv > trimmean(maxenv,percent_to_drop)));
    maxenv(maxenv < residmean) = 0;
    residmaxsmo = conv(maxenv, filttraces,'same');
    
    
    % find the first and second derivatives of the max
    bmax_1st_deriv = diff(residmaxsmo);
    bmax_2nd_deriv = diff(residmaxsmo,2);
    
    % apply a signum filter to get samples at zero crossings and make only 1
    % and 0's
    bsign_1_deriv = sign(bmax_1st_deriv);
    bsign_1_deriv(bsign_1_deriv == 0) = 1;
    
    % find the point where sign of 1st deriv changes
    bdiffsign = diff(bsign_1_deriv);
    bmdiffsign = bdiffsign;
    
    % set the point to zero where second derivative is positive and pad, this
    % finds the peaks in the max dataset
    bdiffsign(sign(bmax_2nd_deriv) > 0) = 0;
    bdiffsign = [1;bdiffsign];
    
    bdiffsign = conv(bdiffsign, passfilt,'same');
    bdiffsign(bdiffsign ~= 0) = 1;
    
    outzone(:,ckk) = bdiffsign;
    
    % % now make the scaler to make envlope all fit the value 2000
    % scalepos = 2000 ./ maxenv;
    % %scalepos = bsxfun(@rdivide,2000,posenv);
    % scalepos(isnan(scalepos)) = 0;
    %
    % % apply a median filter to remove sudden jumps in scaling and the small
    % % averaging to make sure it is a smooth scalar
    % scalepos = medfilt3nt(scalepos,15,0);
    % scalepos =  conv(scalepos,filt_smo,'same');
    %
    %
    % %Apply the scaling to the input data
    % %trim_data_filt(ckk) = bsxfun(@times,scalepos,trim_data(ckk));
    % trim_data_filt(:,ckk) = trim_data(:,ckk).*scalepos;
    
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

