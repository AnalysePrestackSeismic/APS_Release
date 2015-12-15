function trim_data_filt = compress_seismic_to_8_int(trim_data)
%% Defination: TIME_BALENCE scales a gather or setion from the envelope of the amplitudes
% Input:
% Output:
% Writes to Disk:
%scalepeak = single(32700);
scalepeak = double(28000);
%scalepeak = single(110);
morecompress = 1;
envlplot = 0;
outplot = 0;
firstpos = single(0);
endpos = single(0);
firstneg = single(0);
endneg = single(0);
nanperwhtnoise  = 0.02;
if morecompress == 1
    nanperwhtnoise  = 0.1;
    perwhtnoise = 0.3;
elseif morecompress == 2
    perwhtnoise = 0.06;
else
    perwhtnoise = 0.01;
end

%to make the envelope tend to the value +/- 2000
%   Detailed explanation goes here
%%
filt_smo =  ones(1,3)/3;
%filttraces = [1 2 2 3 3 3 2 2 1]/19;
%trim_data_filt = trim_data;
intrlen = size(trim_data,1);
trim_data_filt = zeros(intrlen,size(trim_data,2),'int8');
diffout = zeros(intrlen,1,'single');

%midpointslogic = logical(diffout);
midpoints = diffout;
itpslocs = single(1:intrlen)';
itpslocstest = itpslocs;
origcount = itpslocs;
tmpupper = itpslocstest;

%for ckk = 1:size(trim_data,2)
for ckk = 7:7   
    % reset some variables
    itpslocs = origcount;
    itpslocstest = origcount;
    %tmpupper = origcount;
    % find the variance of the input data
    %varmeas = var(trim_data(:,ckk));
    %fprintf('%-10.8f variance\n',varmeas);
    %trim_data(:,ckk) =  conv(trim_data(:,ckk),filttraces,'same');
    if envlplot == 1
        figure(2); plot(trim_data(:,ckk));
        hold on
    end
    % find the first and second derivatives of the max
    max_1st_deriv = diff(trim_data(:,ckk));
    max_2nd_deriv = diff(trim_data(:,ckk),2);
    
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
    itpsval = trim_data((diffsign < 0),ckk);
    %itpsval = itpsval + (trimmean(itpsval,10)*0.001);
    itpslocsin = itpslocs(diffsign < 0);
    
    mitpsval = trim_data((mdiffsign > 0),ckk);
    %mitpslocs = single(1:size(trim_data,1))';
    mitpslocsin = itpslocs(mdiffsign > 0);
    
    if itpslocsin(end) > mitpslocsin(end)
        if itpslocsin(end) ~= intrlen
            itpslocsin = [itpslocsin; intrlen];
            itpsval = [itpsval; itpsval(end)];
        end
        mitpslocsin = [mitpslocsin; intrlen];
        mitpsval = [mitpsval; trim_data(end,ckk)];
    else
        if mitpslocsin(end) ~= intrlen
            mitpslocsin = [mitpslocsin; intrlen];
            mitpsval = [mitpsval; mitpsval(end)];
        end
        itpslocsin = [itpslocsin; intrlen];
        itpsval = [itpsval; trim_data(end,ckk)];
    end
    
    if itpslocsin(1) < mitpslocsin(1)
        if itpslocsin(1) ~= 1
            itpslocsin = [1; itpslocsin];
            itpsval = [itpsval(1); itpsval];
        end
        mitpslocsin = [1; mitpslocsin];
        mitpsval = [trim_data(1,ckk); mitpsval];
    else
        if mitpslocsin(1) ~= 1
            mitpslocsin = [1; mitpslocsin];
            mitpsval = [mitpsval(1); mitpsval];
        end
        itpslocsin = [1; itpslocsin];
        itpsval = [trim_data(1,ckk); itpsval];
    end
    
    %store orginal points to compare at the end
    origuprloc = itpslocsin;
    origuprval = itpsval;
    origlowloc = mitpslocsin;
    origlowval = mitpsval;    

    % make sure the start and end do not get lost
    firstpos = itpsval(1);
    endpos = itpsval(end);
    firstneg = mitpsval(1);
    endneg = mitpsval(end);
 
%     if itpslocsin(end) ~= intrlen
%         itpslocsin = [itpslocsin; intrlen];
%         itpsval = [itpsval; endpos];
%     end
%     if itpslocsin(1) ~= 1
%         itpslocsin = [1; itpslocsin];
%         itpsval = [firstpos; itpsval];
%     end
% 
%     %mitpslocsin(end) = intrlen;
%     %mitpslocsin(1) = 1;
% 
%     if mitpslocsin(end) ~= intrlen
%         mitpslocsin = [mitpslocsin; intrlen];
%         mitpsval = [mitpsval; endneg];
%     end
%     if mitpslocsin(1) ~= 1
%         mitpslocsin = [1; mitpslocsin];
%         mitpsval = [firstneg; mitpsval];
%     end    
       
    %itpslocsin(end) = intrlen;
    %itpslocsin(1) = 1;
    % interpolate to make the envelope only using fast linear interp
    posenv = double(interp1q(itpslocsin,itpsval,origcount));
    if envlplot == 1
        plot(posenv,'-r')
    end
    if (morecompress > 1)
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
        % make sure the start and end do not get lost
        if itpslocsin(end) ~= intrlen
            itpslocsin = [itpslocsin; intrlen];
            itpsval = [itpsval; endpos];
        end
        if itpslocsin(1) ~= 1
            itpslocsin = [1; itpslocsin];
            itpsval = [firstpos; itpsval];
        end
        
        % add small white noise to avoid division errors/mess
        %itpsval = itpsval + (trimmean(itpsval,10)*0.001);
        % interpolate to make the envelope only using fast linear interp
        posenv = double(interp1q(itpslocsin,itpsval,origcount));
        %==========================================================================
        %test to see if any points have been chopped off in the envelope
        newposvals =  posenv(origuprloc);
        tmpupper = origcount;
        tmpupper(itpslocsin) = itpsval;
        tmpupper(origuprloc(newposvals < origuprval)) = origuprval(newposvals < origuprval);
        itpslocsin = itpslocstest(tmpupper ~= itpslocstest);
        itpsval = tmpupper(tmpupper ~= itpslocstest);
        posenv = double(interp1q(itpslocsin,itpsval,origcount));
        %==========================================================================
        if envlplot == 1
            plot(posenv,'-g')
        end
        %==========================================================================
    end
    if (morecompress > 2)
        %==========================================================================
        % calc the envelope of the envelope of the envelope
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
        
        % make sure the start and end do not get lost
        if itpslocsin(end) ~= intrlen
            itpslocsin = [itpslocsin; intrlen];
            itpsval = [itpsval; endpos];
        end
        if itpslocsin(1) ~= 1
            itpslocsin = [1; itpslocsin];
            itpsval = [firstpos; itpsval];
        end
        
        posenv = double(interp1q(itpslocsin,itpsval,origcount));  
        %==========================================================================
        %test to see if any points have been chopped off in the envelope
        newposvals =  posenv(origuprloc);
        tmpupper = origcount;
        tmpupper(itpslocsin) = itpsval;
        tmpupper(origuprloc(newposvals < origuprval)) = origuprval(newposvals < origuprval);
        itpslocsin = itpslocstest(tmpupper ~= itpslocstest);
        itpsval = tmpupper(tmpupper ~= itpslocstest);
        
        %==========================================================================
        
        
        %==========================================================================
    end
    % add start and end values
    %itpslocsin(end) = intrlen;
    %itpslocsin(1) = 1;
    % add small white noise to avoid division errors/mess
    itpsval = itpsval + (abs(trimmean(itpsval,20))*nanperwhtnoise);
    % interpolate to make the envelope only using fast linear interp
    posenv = double(interp1q(itpslocsin,itpsval,origcount));       
    if envlplot == 1      
        plot(posenv,'-m')
    end
    
    
    %==========================================================================
    %use the mins logical to get the values and indexes to then interpolate to
    %make and envelope which includes the signal, but preserves the wiggles in
    %the dataset
%     mitpsval = trim_data((orig_mdiffsign > 0),ckk);
%     %mitpslocs = single(1:size(trim_data,1))';
%     mitpslocsin = itpslocs(orig_mdiffsign > 0);
%     
%     %mitpslocsin(end) = intrlen;
%     %mitpslocsin(1) = 1;
%     firstneg = mitpsval(1);
%     endneg = mitpsval(end);
%     if mitpslocsin(end) ~= intrlen
%         mitpslocsin = [mitpslocsin; intrlen];
%         mitpsval = [mitpsval; endneg];
%     end
%     if mitpslocsin(1) ~= 1
%         mitpslocsin = [1; mitpslocsin];
%         mitpsval = [firstneg; mitpsval];
%     end    
    
    %mitpsval = mitpsval + (trimmean(mitpsval,10)*0.001);
    % interpolate to make the min envelope only using fast linear interp
    mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
    if envlplot == 1
        plot(mposenv,'-r')
    end
    if (morecompress > 1)
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
        mitpsval = mposenv(mdiffsign > 0);
        mitpslocsin = itpslocs(mdiffsign > 0);
        % add small white noise to avoid division errors/mess
        %mitpsval = mitpsval + (trimmean(mitpsval,10)*0.001);
        % interpolate to make the envelope only using fast linear interp
        
        % make sure the start and end do not get lost
        if mitpslocsin(end) ~= intrlen
            mitpslocsin = [mitpslocsin; intrlen];
            mitpsval = [mitpsval; endneg];
        end
        if mitpslocsin(1) ~= 1
            mitpslocsin = [1; mitpslocsin];
            mitpsval = [firstneg; mitpsval];
        end
        
        mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
        
        %==========================================================================
        %test to see if any points have been chopped off in the envelope
        newposvals =  mposenv(origlowloc);
        tmpupper = origcount;
        tmpupper(mitpslocsin) = mitpsval;
        tmpupper(origlowloc(newposvals > origlowval)) = origlowval(newposvals > origlowval);
        mitpslocsin = itpslocstest(tmpupper ~= itpslocstest);
        mitpsval = tmpupper(tmpupper ~= itpslocstest);
        mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
        %==========================================================================        
        
        if envlplot == 1
            plot(mposenv,'-g')
        end
    end
    %==========================================================================
    if (morecompress > 2)
        % calc the envelope of the envelope of the envelope
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
        mitpsval = mposenv(mdiffsign > 0);
        mitpslocsin = itpslocs(mdiffsign > 0);
        % make sure the start and end do not get lost
        if mitpslocsin(end) ~= intrlen
            mitpslocsin = [mitpslocsin; intrlen];
            mitpsval = [mitpsval; endneg];
        end
        if mitpslocsin(1) ~= 1
            mitpslocsin = [1; mitpslocsin];
            mitpsval = [firstneg; mitpsval];
        end
        mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
        %==========================================================================
        %test to see if any points have been chopped off in the envelope
        newposvals =  mposenv(origlowloc);
        tmpupper = origcount;
        tmpupper(mitpslocsin) = mitpsval;
        tmpupper(origlowloc(newposvals > origlowval)) = origlowval(newposvals > origlowval);
        mitpslocsin = itpslocstest(tmpupper ~= itpslocstest);
        mitpsval = tmpupper(tmpupper ~= itpslocstest);
        mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));      
        
        %==========================================================================
    end
    % add small white noise to avoid division errors/mess
    % add start and end values
    %mitpslocsin(end) = intrlen;
    %mitpslocsin(1) = 1;    
    mitpsval = mitpsval - (abs(trimmean(mitpsval,20))*nanperwhtnoise);
    % interpolate to make the envelope only using fast linear interp
    mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
    
    if envlplot == 1
        plot(mposenv-(abs((posenv - mposenv)).*perwhtnoise),'-om')
        plot(posenv+(abs((posenv - mposenv)).*perwhtnoise),'-om')
    end
    %=========================================================================

    % add % white noise based on size of wiggle
    posenv = posenv+(abs((posenv - mposenv)).*perwhtnoise);
    mposenv = mposenv-(abs((posenv - mposenv)).*perwhtnoise);
    
    midpoints = (posenv - mposenv).*0.5;
    scalepos = scalepeak ./ midpoints;
    %itpsval = abs(single(scalepeak ./ itpsval));
    %mitpsval = abs(single(-scalepeak ./ mitpsval));
    
    
    %itpslocs(mitpslocsin) = mitpsval;
    %itpslocs(itpslocsin) = itpsval;
    %itpslocs(1) = scalepeak;    % add start and end values to avoid missing points
    %itpslocs(end) = itpslocs(1);
    %alllocsin = itpslocstest(itpslocs ~= itpslocstest);
    %allvalsin = itpslocs(itpslocs ~= itpslocstest);
    
    %allvalsin(isnan(allvalsin)) = 0;
    
    %mitpsval(end+1) = scalepeak;    % add start and end values to avoid missing points
    %mitpslocsin(end) = intrlen;
    %itpsval(end+1) = scalepeak;    % add start and end values to avoid missing points
    %itpslocsin(end) = intrlen;
    
    if envlplot == 1
        plot(itpslocsin,(itpsval),'--oc');
        plot(mitpslocsin,(mitpsval),'--og');
    end
    % interpolate to make the min envelope only using fast linear interp
    %scalepos = single(interp1q(alllocsin,allvalsin,itpslocstest));
    
    %scalepos = single(interp1q(itpslocsin,itpsval,itpslocstest));
    %scaleneg = single(interp1q(mitpslocsin,mitpsval,itpslocstest));
    
    if envlplot == 1
        figure(3); plot(scalepos);
        hold on;
        %plot(scaleneg,'-r');
    end
    %maxenv = max([abs(posenv) abs(mposenv)],[],2);
    
    
    % now make the scaler to make envlope all fit the value 2000
    %scalepos = 2000 ./ maxenv;
    %scalepos = bsxfun(@rdivide,2000,posenv);
    scalepos(isnan(scalepos)) = 1;
    scalepos_restore = single(1./scalepos);
    scalepos_restore(isnan(scalepos_restore)) = 1;
    
    %scaleneg(isnan(scaleneg)) = -1;
    %%scaleneg = scaleneg*-1;
    %scaleneg_restore = single(1./scaleneg);
    %scaleneg_restore(isnan(scaleneg_restore)) = 1;
    
    scalemean = trimmean(scalepos,50);
    
    % apply a median filter to remove sudden jumps in scaling and the small
    % averaging to make sure it is a smooth scalar
    %scalepos = medfilt3nt(scalepos,15,0);
    %scalepos =  conv(scalepos,filt_smo,'same');
    
    
    %Apply the scaling to the input data
    %%trim_data_filt(ckk) = bsxfun(@times,scalepos,trim_data(ckk));
    
    trim_data_filt(:,ckk) = int8((trim_data(:,ckk)-(mposenv+midpoints)).*scalepos);
    
    %midpointslogic = (trim_data(:,ckk) > (mposenv + (posenv - mposenv)/2));
    %midpointslogic = (trim_data(:,ckk) > 0);

    
    %trim_data_filt(midpointslogic,ckk) = int16(trim_data(midpointslogic,ckk).*scalepos(midpointslogic));
    %trim_data_filt((midpointslogic == 0),ckk) = int16(trim_data((midpointslogic == 0),ckk).*scaleneg(midpointslogic == 0));
    %%%trim_data_filt((trim_data(:,ckk) < 0),ckk) = int16(trim_data((trim_data(:,ckk) < 0),ckk).*scaleneg(trim_data(:,ckk) < 0));

    
    if envlplot == 1
        figure(36); plot(posenv)
        hold on;
        plot(mposenv,'-r')
        plot((mposenv+midpoints),'-g')
        plot(trim_data_filt(:,ckk),'-c')
        plot(trim_data(:,ckk),'-m');
        %plotlogic = (((midpointslogic *2)+1)-2)*1000;
        %plot(plotlogic,'--g')
        plot(scalepos.*500,'--r')
    end    
    %figure(4); plot(trim_data_filt(:,ckk));
    
    
    
    %=====================================================================
    %uncompress the data to compare
    
    diffout = single((single(trim_data_filt(:,ckk)).*scalepos_restore)+(mposenv+midpoints));
    
    %diffout(trim_data_filt(:,ckk)>0) = single(single(trim_data_filt((trim_data_filt(:,ckk)>0),ckk)).*scalepos_restore(trim_data_filt(:,ckk)>0));
    %diffout(trim_data_filt(:,ckk)<0) = single(single(trim_data_filt((trim_data_filt(:,ckk)<0),ckk)).*scaleneg_restore(trim_data_filt(:,ckk)<0));
    
    %diffout(midpointslogic) = single(single(trim_data_filt(midpointslogic,ckk)).*scalepos_restore(midpointslogic));
    %diffout(midpointslogic == 0) = single(single(trim_data_filt((midpointslogic == 0),ckk)).*scaleneg_restore(midpointslogic == 0));

    
    uncerror =  (trim_data(:,ckk) - diffout) ./ trim_data(:,ckk);
    uncerror(isnan(uncerror)) = 0;
    uncerror(isinf(uncerror)) = 0;
    
    summoferror = sum(abs((trim_data(:,ckk) - diffout)));
    summofinput = sum(abs(trim_data(:,ckk)));
    fprintf('%-10.8f percent error in reconstruction of dataset, one in %d, using %d scalars\n',(summoferror/summofinput)*100,round((summofinput/summoferror)),(size(mitpsval,1)*4));

    uncerror(scalepos > scalemean*5) = 0;
    unerrcount = sum(uncerror ~= 0);
    avgpercterr = (sum(abs(uncerror)))/unerrcount;

    %fprintf('%-10.8f percent error in reconstruction of dataset, one in %d, using %d scalars\n',avgpercterr,round(100/avgpercterr),(size(mitpsval,1)*4));
    if outplot == 1
        figure; plot(diffout);
        %hold on
        figure; plot(trim_data(:,ckk),'-g');
    end
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

