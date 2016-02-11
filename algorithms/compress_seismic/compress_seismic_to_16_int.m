function trim_data_filt = compress_seismic_to_16_int( trim_data )
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.
%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
%% Defination: TIME_BALENCE scales a gather or setion from the envelope of the amplitudes
% Input:
% Output:
% Writes to Disk:
%scalepeak = single(32700);
scalepeak = double(28000);
%scalepeak = single(110);
morecompress = 3;
envlplot = 4;
outplot = 0;
firstpos = single(0);
endpos = single(0);
firstneg = single(0);
endneg = single(0);
%nanperwhtnoise  = 0.02;
nanperwhtnoise  = 1.0001;
if morecompress == 1
    nanperwhtnoise  = 1.0001;
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
trim_data_filt = zeros(intrlen,size(trim_data,2),'int16');
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
    
    %itpsval = itpsval + abs(mean(itpsval)*nanperwhtnoise);
%    itpsval = itpsval + (abs(trimmean(itpsval,20))*nanperwhtnoise);
    itpsval = itpsval.*nanperwhtnoise;
    
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
        %mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));      
        
        %==========================================================================
    end
    % add small white noise to avoid division errors/mess
    % add start and end values
    %mitpslocsin(end) = intrlen;
    %mitpslocsin(1) = 1;    
    
    %was adding some white noise
    %mitpsval = mitpsval - (abs(mean(mitpsval))*nanperwhtnoise);
%    mitpsval = mitpsval - (abs(trimmean(mitpsval,20))*nanperwhtnoise);
    mitpsval = mitpsval.*nanperwhtnoise;
    
    % interpolate to make the envelope only using fast linear interp
    mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
    
    if envlplot == 1
        plot(mposenv-(abs((posenv - mposenv)).*perwhtnoise),'-om')
        plot(posenv+(abs((posenv - mposenv)).*perwhtnoise),'-om')
    end
    %=========================================================================

    % add % white noise based on size of wiggle
    %posenv = posenv+(abs((posenv - mposenv)).*perwhtnoise);
    %mposenv = mposenv-(abs((posenv - mposenv)).*perwhtnoise);
    
    % calculate the mid points of the upper and lower envelopes from the
    % interpolated envelopes
    midpoints = posenv+(((posenv.*-1) - (mposenv.*-1)).*0.5);
    
    %decimate the mid points and then reinterpolate same as the
    %uncompression would do
    midpoints_red = midpoints(mitpslocsin);
    midpoints = double(interp1q(mitpslocsin,midpoints_red,origcount));
    
    scalepos = scalepeak ./ midpoints;
    
    scalepos(isnan(scalepos)) = 1;
    scalepos_restore = single(1./scalepos);
    scalepos_restore(isnan(scalepos_restore)) = 1;
    
    if envlplot == 2
        figure(33); plot(posenv,'-k');
        hold on;
        plot(mposenv,'-r');
        plot(midpoints,'-b');
        plot(mitpslocsin,midpoints_red,'-g')
        figure(34); plot(scalepos);
    end
    
    
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

    
    %scaleneg(isnan(scaleneg)) = -1;
    %%scaleneg = scaleneg*-1;
    %scaleneg_restore = single(1./scaleneg);
    %scaleneg_restore(isnan(scaleneg_restore)) = 1;
    
%    scalemean = trimmean(scalepos,50);
    
    % apply a median filter to remove sudden jumps in scaling and the small
    % averaging to make sure it is a smooth scalar
    %scalepos = medfilt3nt(scalepos,15,0);
    %scalepos =  conv(scalepos,filt_smo,'same');
    
    
    
    if envlplot == 3
       
        figure(97); plot(trim_data(:,ckk),'-k');
        hold on;
        plot(mposenv,'-r');
        plot((midpoints-(mposenv-midpoints)),'-g');
        plot(midpoints,'-b');
        %plot((midpoints - mposenv),'-m');
    end
    
    
    %Apply the scaling to the input data
    %%trim_data_filt(ckk) = bsxfun(@times,scalepos,trim_data(ckk));
    
    %trim_data_filt(:,ckk) = int16((trim_data(:,ckk)-(mposenv+midpoints)).*scalepos);
    scaleval = scalepeak./(midpoints - mposenv);
    scaleval(isnan(scaleval)) = 1;
    trim_data_filt(:,ckk) = int16((trim_data(:,ckk)-midpoints).*scaleval);
    %midpointslogic = (trim_data(:,ckk) > (mposenv + (posenv - mposenv)/2));
    %midpointslogic = (trim_data(:,ckk) > 0);

     if envlplot == 3
         figure(98); plot(trim_data_filt(:,ckk),'-k');
     end
         
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
        %plot(scalepos.*500,'--r')
    end    
    %figure(4); plot(trim_data_filt(:,ckk));
    
   
    
    %=====================================================================
    %uncompress the data to compare
    
    % need to have saved 
%     1) midpoints_red
%     2) mitpslocsin
%     3) mitpsval
%     4) scalepeak
    
    
    samples_arr = single(1:size(trim_data,1))';
    resmposenv = double(interp1q(mitpslocsin,mitpsval,samples_arr));

    res_midpoints = double(interp1q(mitpslocsin,midpoints_red,samples_arr));

    scaleval_res = (res_midpoints - resmposenv)./scalepeak;
    scaleval_res(isnan(scaleval_res)) = 1;

    diffout = single((single(trim_data_filt(:,ckk)).*scaleval_res)+res_midpoints);

    
%     uncerror =  (trim_data(:,ckk) - diffout) ./ trim_data(:,ckk);
%     uncerror(isnan(uncerror)) = 0;
%     uncerror(isinf(uncerror)) = 0;
    
    summoferror = sum(abs((trim_data(:,ckk) - diffout)));
    summofinput = sum(abs(trim_data(:,ckk)));
    fprintf('%-10.8f percent error in reconstruction of dataset, one in %d, using %d scalars\n',(summoferror/summofinput)*100,round((summofinput/summoferror)),((size(mitpsval,1)*3)+1));

%    uncerror(scalepos > scalemean*5) = 0;
%    unerrcount = sum(uncerror ~= 0);
%    avgpercterr = (sum(abs(uncerror)))/unerrcount;

    %fprintf('%-10.8f percent error in reconstruction of dataset, one in %d, using %d scalars\n',avgpercterr,round(100/avgpercterr),(size(mitpsval,1)*4));
    if outplot == 1
        figure(50); plot(diffout);
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

