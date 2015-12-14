function [ trim_data_filt, scalars_out, zlocs_out, tpgrad_orig, scalepeak,  orig_nsamples ] = compress_seis_to_16bit_int( trim_data, morecompress, envlplot )
%
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
% ------------------ License  ------------------
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
% github
% https://github.com/AnalysePrestackSeismic/
% ------------------ FUNCTION DEFINITION ---------------------------------
%
%
% Defination: compress_seis_to_16_or_8_bit_int compress the seismic into 16 bit
% integer or 8 bit with scalars to provide good recovery whilst still looking like
% seismic, so the result cane be viewed at 16 bit int, 8 bit int(like a very short agc)
% or uncompress to 32 bit ieee float
% morecompress values
% 1 = first envelope for data velaue  - do not use very inefficent compression
% 2 = second envelope, envelope of first envelope - 8 bit integer ~ 3 times dependant on trace length
% 3 = third envelope, ie envelope of 2nd envelope - 16 bit integer ~ 1.9 times compression
%
% Input: 
% trim_data = a 2 dimensional array, rows are samples , columns are traces
% morecompress = 2 = 8 bit, 3 = 16 bit output
% envlplot = only set to 3 will produce plots all other plots commented out
%
% Output: 
% trim_data_filt = the output 2 dimensional array of int16 or int 8 numbers
% scalars_out = 3 d array of the scalars, size of 2; number of samples; number of traces
%       1st dimension is the mid points (e.g. scalars_out(1,:,:) , second dimension is the
%       lower envelope e.g. scalars_out(2,:,:), values are singles
% zlocs_out = the z locations out 2d array, (locations, traces) stored as
%       unit16
% tpgrad_orig = number of scalar points
% scalepeak = peakvalue scaled to
% orig_nsamples = orginal number of samples in the trace
%
%
% to save space the code only write the sample numnbers as uint16 so means
% that the trace length limit is 65536 samples, which is 131 seconds at 2ms
% sampling, could just change to uint32 to move to longer records, but that
% would also need a review of the number of scalars
%
% Authors: Charles Jones 2015
%%
%
% read the command line variables
morecompress = str2double(morecompress);
envlplot = str2double(envlplot);
% work out input dimensions
no_of_traces = size(trim_data,2);
intrlen = size(trim_data,1);
orig_nsamples = uint16(intrlen);
%
outplot = 0;
firstpos = single(0);
endpos = single(0);
firstneg = single(0);
endneg = single(0);
%nanperwhtnoise  = 0.02;

if morecompress == 1
    tpgrad_orig = floor(intrlen/11);  % max number of scalars to store
    nanperwhtnoise  = 1.15;  % % to add to avoid cliping
    %perwhtnoise = 0.3;
    scalepeak = double(120);   % the max value to scale to
    smalldrop = 0.0000005;
elseif morecompress == 2
    tpgrad_orig = floor(intrlen/33);   % max number of scalars to store
    nanperwhtnoise  = 1.05;  % % to add to avoid cliping
    %perwhtnoise = 0.06;
    scalepeak = double(120); % the max value to scale to
    smalldrop = 0.00005;
else
    tpgrad_orig = floor(intrlen/66);   % max number of scalars to store
    nanperwhtnoise  = 1.05;  % % to add to avoid cliping
    %perwhtnoise = 0.01;
    scalepeak = double(28000); % the max value to scale to
    smalldrop = 0.00001;
end

%to make the envelope tend to the value +/- 2000
%   Detailed explanation goes here
%%
midpoints_red_t = [ones(tpgrad_orig,1,'single') zeros(tpgrad_orig,1,'single')];
midpoints_red_torg = midpoints_red_t;
%
%filt_smo =  ones(1,3)/3;
%filttraces = [1 2 2 3 3 3 2 2 1]/19;
%trim_data_filt = trim_data;


if morecompress < 3
    trim_data_filt = zeros(intrlen,no_of_traces,'int8');
else
    trim_data_filt = zeros(intrlen,no_of_traces,'int16');
end
diffout = zeros(intrlen,1,'single');
midblank = zeros(intrlen,1,'double');
scaleval = midblank;
%midpointslogic = logical(diffout);
midpoints = midblank;
itpslocs = single(1:intrlen)';
itpslocstest = itpslocs;
origcount = itpslocs;
tmpupper = itpslocstest;
blankone = single(1);
blanklogic = logical(1);
blankdbl = 0;
%
scalars_out = ones(2,tpgrad_orig,no_of_traces,'single');
zlocs_out = ones(tpgrad_orig,no_of_traces,'uint16');
%
%
%
for ckk = 1:no_of_traces
    %for ckk = 1227:1227
    % reset some variables
    itpslocs = origcount;
    itpslocstest = origcount;
    midpoints_red_t = midpoints_red_torg;
    itpslocsin = blankone;
    itpsval = blankone;
    mid_lens = blankone;
    midpoints_red  = blankone;
    midpoints_red_grad  = blankone;
    mitpslocsin  = blankone;
    mitpsval = blankone;
    mitpsval_red = blankone;
    points_all = blankone;
    top15grads = blankone;
    xgrad = blankone;
    ygrad = blankone;
    
    summoferror = blankdbl;
    summofinput = blankdbl;
    pointlogic = blanklogic;
    tpgrad = tpgrad_orig;
    %tmpupper = origcount;
    % find the variance of the input data
    %varmeas = var(trim_data(:,ckk));
    %fprintf('%-10.8f variance\n',varmeas);
    %trim_data(:,ckk) =  conv(trim_data(:,ckk),filttraces,'same');
    %     if envlplot == 1
    %         figure(2); plot(trim_data(:,ckk));
    %         hold on
    %     end
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
    
    % interpolate to make the envelope only using fast linear interp
    %posenv = double(interp1q(itpslocsin,itpsval,origcount));
    posenv = double(makefastinterp1(double(itpslocsin),double(itpsval),double(origcount)));
    
    %     if envlplot == 1
    %         plot(posenv,'-r')
    %     end
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
        %posenv = double(interp1q(itpslocsin,itpsval,origcount));
        posenv = double(makefastinterp1(double(itpslocsin),double(itpsval),double(origcount)));
        %==========================================================================
        %test to see if any points have been chopped off in the envelope
        newposvals =  posenv(origuprloc);
        tmpupper = origcount;
        tmpupper(itpslocsin) = itpsval;
        tmpupper(origuprloc(newposvals < origuprval)) = origuprval(newposvals < origuprval);
        itpslocsin = itpslocstest(tmpupper ~= itpslocstest);
        itpsval = tmpupper(tmpupper ~= itpslocstest);
        %posenv = double(interp1q(itpslocsin,itpsval,origcount));
        posenv = double(makefastinterp1(double(itpslocsin),double(itpsval),double(origcount)));
        %==========================================================================
        %         if envlplot == 1
        %             plot(posenv,'-g')
        %         end
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
        
        %posenv = double(interp1q(itpslocsin,itpsval,origcount));
        posenv = double(makefastinterp1(double(itpslocsin),double(itpsval),double(origcount)));
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
    
    finallocs_upper = ((itpsval.*(abs(itpsval -[0;itpsval(1:(end-1))]) + abs(itpsval -[itpsval(2:(end));0])))./itpsval) > 0.00001;
    finallocs_upper(1:2) = 1;
    finallocs_upper(end) = 1;
    itpslocsin = itpslocsin(finallocs_upper);
    itpsval = itpsval(finallocs_upper);
    
    
    % add small white noise to avoid division errors/mess
    itpsval = itpsval.*nanperwhtnoise;
    
    % interpolate to make the envelope only using fast linear interp
    %posenv = double(interp1q(itpslocsin,itpsval,origcount));
    posenv = double(makefastinterp1(double(itpslocsin),double(itpsval),double(origcount)));
    %mposenv = double(makefastinterp1(double(mitpslocsin),double(mitpsval),double(origcount)));
    %     if envlplot == 1
    %         plot(posenv,'-m')
    %     end
    
    
    %==========================================================================
    %use the mins logical to get the values and indexes to then interpolate to
    %make and envelope which includes the signal, but preserves the wiggles in
    %the dataset
    
    % interpolate to make the min envelope only using fast linear interp
    %mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
    mposenv = double(makefastinterp1(double(mitpslocsin),double(mitpsval),double(origcount)));
    %     if envlplot == 1
    %         plot(mposenv,'-r')
    %     end
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
        
        %mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
        mposenv = double(makefastinterp1(double(mitpslocsin),double(mitpsval),double(origcount)));
        %==========================================================================
        %test to see if any points have been chopped off in the envelope
        newposvals =  mposenv(origlowloc);
        tmpupper = origcount;
        tmpupper(mitpslocsin) = mitpsval;
        tmpupper(origlowloc(newposvals > origlowval)) = origlowval(newposvals > origlowval);
        mitpslocsin = itpslocstest(tmpupper ~= itpslocstest);
        mitpsval = tmpupper(tmpupper ~= itpslocstest);
        %mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
        mposenv = double(makefastinterp1(double(mitpslocsin),double(mitpsval),double(origcount)));
        %==========================================================================
        
        %         if envlplot == 1
        %             plot(mposenv,'-g')
        %         end
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
        %mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
        mposenv = double(makefastinterp1(double(mitpslocsin),double(mitpsval),double(origcount)));
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
    % remove points that are very nearly the same value
    %finallocs = [1;mitpslocsin(((mitpsval.*(abs(mitpsval -[0;mitpsval(1:(end-1))]) + abs(mitpsval -[mitpsval(2:(end));0])))./mitpsval) > 0.00001)];
    finallocs_lower = ((mitpsval.*(abs(mitpsval -[0;mitpsval(1:(end-1))]) + abs(mitpsval -[mitpsval(2:(end));0])))./mitpsval) > 0.00001;
    finallocs_lower(1:2) = 1;
    finallocs_lower(end) = 1;
    mitpslocsin = mitpslocsin(finallocs_lower);
    mitpsval = mitpsval(finallocs_lower);
    
    %was adding some white noise
    %    mitpsval = mitpsval - (abs(trimmean(mitpsval,20))*nanperwhtnoise);
    mitpsval = mitpsval.*nanperwhtnoise;
    
    % interpolate to make the envelope only using fast linear interp
    %mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
    mposenv = double(makefastinterp1(double(mitpslocsin),double(mitpsval),double(origcount)));
    
    
    %     if envlplot == 1
    %         plot(mposenv-(abs((posenv - mposenv)).*perwhtnoise),'-om')
    %         plot(posenv+(abs((posenv - mposenv)).*perwhtnoise),'-om')
    %     end
    %=========================================================================
    
    % add % white noise based on size of wiggle
    %posenv = posenv+(abs((posenv - mposenv)).*perwhtnoise);
    %mposenv = mposenv-(abs((posenv - mposenv)).*perwhtnoise);
    
    % calculate the mid points of the upper and lower envelopes from the
    % interpolated envelopes
    midpoints = posenv+(((posenv.*-1) - (mposenv.*-1)).*0.5);
    
    %decimate the mid points and then reinterpolate same as the
    %uncompression would do
    %if morecompress < 3
    % make the mid points for all the points on up and lower envelope
    itpslocs(itpslocsin) = itpsval;
    itpslocs(mitpslocsin) = mitpsval;
    points_all = itpslocstest(itpslocs ~= itpslocstest);
    %else
    %    points_all = mitpslocsin;
    %end
    
    
    % make the sets of points to store
    mitpsval_red = single(mposenv(points_all));
    midpoints_red = single(midpoints(points_all));
    
    %drop points that are next to each other
    ygrad = ([midpoints_red(2:end);midpoints_red(end)] - midpoints_red);
    xgrad = ( [points_all(2:end);points_all(end)]  - points_all   );
    
    
    if envlplot == 3
        figure(22); plot(points_all,midpoints_red,'-ob'); hold on
        plot(points_all,mitpsval_red,'-ob')
        plot((midpoints-(mposenv-midpoints)),'-g');
        plot(trim_data(:,ckk),'-k')
    end
    
    
    % ====================================================================
    % decimate mid points to just those points which are not too close,
    
    % but keep start and end
    %       points_all = points_all([999;999;diff(points_all(2:end-1));999]>4);
    % workout the distance between the points and replace those too close to
    % each other with the average, ie on do test for those point less than 10 samples apart
    mid_lens = (ygrad.*ygrad) + (xgrad.*xgrad);
    mean_mid_lens = mean(mid_lens);
    pointlogic =or([999;999;diff(points_all(3:end-1));999;999]>6, mid_lens > (mean_mid_lens*smalldrop));
    points_all = points_all(pointlogic);
    ygrad = ygrad(pointlogic);
    xgrad = xgrad(pointlogic);
    
    % re-make the sets of points to store
    mitpsval_red = single(mposenv(points_all));
    midpoints_red = single(midpoints(points_all));
    
    % so instread calculate grdients of the mid points and selct those with
    % similar gradients either side
    %midpoints_red_grad = ([midpoints_red(2:end);midpoints_red(end)] - midpoints_red)./( [points_all(2:end);points_all(end)]  - points_all   );
    midpoints_red_grad = ygrad./xgrad;
    
    %bob = [midpoints_red [1;midpoints_red_grad(1:(end-1))] ([1;midpoints_red_grad(1:(end-1))]-[midpoints_red_grad(1:(end-1));1])  abs([1;midpoints_red_grad(1:(end-1))]-[midpoints_red_grad(1:(end-1));1])  ];
    top15grads = sortrows([abs([9999999;midpoints_red_grad(1:(end-2));9999999]-[midpoints_red_grad(1:(end-2));9999999;0])  points_all ],-1);
    
    %populate the fixed size array to hold the top n points
    if tpgrad > size(top15grads,1);
        tpgrad = size(top15grads,1);
    end
    midpoints_red_t(1:tpgrad,:) = [top15grads(1:tpgrad,2) single(midpoints(top15grads(1:tpgrad,2)))];
    midpoints_red_t = sortrows(midpoints_red_t);
    
    
    if envlplot == 3
        figure(1); plot(points_all,midpoints_red,'-ob'); hold on
        plot(midpoints_red_t(:,1),midpoints_red_t(:,2),'or');
        plot(points_all,mitpsval_red,'-ob')
        plot(midpoints_red_t(:,1),mposenv(midpoints_red_t(:,1)),'-or');
        plot((midpoints-(mposenv-midpoints)),'-g');
        plot(trim_data(:,ckk),'-k')
    end
    
    midpoints_red = midpoints_red_t(:,2);
    points_all = midpoints_red_t(:,1);
    mitpsval_red = single(mposenv(points_all));
    
    
    %     %test to see if any points have been chopped off in the envelope
    %     newposvals =  mposenv(origlowloc);
    %     tmpupper = origcount;
    %     tmpupper(mitpslocsin) = mitpsval;
    %     tmpupper(origlowloc(newposvals > origlowval)) = origlowval(newposvals > origlowval);
    %     mitpslocsin = itpslocstest(tmpupper ~= itpslocstest);
    %     mitpsval = tmpupper(tmpupper ~= itpslocstest);
    %     mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
    %
    %
    %     % remove points that are very nearly the same value
    %     finallocs_lower = ((midpoints_red.*(abs(midpoints_red -[0;midpoints_red(1:(end-1))]) + abs(midpoints_red -[midpoints_red(2:(end));0])))./midpoints_red) > 0.001;
    %     finallocs_lower(1:2) = 1;
    %     finallocs_lower(end) = 1;
    %     points_all = points_all(finallocs_lower);
    %     midpoints_red = midpoints_red(finallocs_lower);
    %     mitpsval_red = mitpsval_red(finallocs_lower);
    
    
    %=====================================================================
    % re interp to be consistent with decompression
    %    midpoints = double(interp1q(points_all,midpoints_red,origcount));
    %    mposenv = double(interp1q(points_all,mitpsval_red,origcount));
    
    % a different version using a mex file for the interp
    midpoints = double(makefastinterp1(double(points_all),double(midpoints_red),double(origcount)));
    mposenv = double(makefastinterp1(double(points_all),double(mitpsval_red),double(origcount)));
    
    % now make the scaler to make envlope all fit the value 2000
    %scalepos = 2000 ./ maxenv;
    %scalepos = bsxfun(@rdivide,2000,posenv);
    
    if envlplot == 3
        figure(97); plot(trim_data(:,ckk),'-k');
        hold on;
        plot(mposenv,'-r');
        plot((midpoints-(mposenv-midpoints)),'-g');
        %plot(posenv,'-m');
        plot(midpoints,'-b');
        %plot((midpoints - mposenv),'-m');
    end
    
    
    %Apply the scaling to the input data
    scaleval = scalepeak./(midpoints - mposenv);
    scaleval(isnan(scaleval)) = 1;
    scaleval(isinf(scaleval)) = 1;
    
    if morecompress < 3
        trim_data_filt(:,ckk) = int8((trim_data(:,ckk)-midpoints).*scaleval);
    else
        trim_data_filt(:,ckk) = int16((trim_data(:,ckk)-midpoints).*scaleval);
    end
    
    if envlplot == 3
        figure(99); plot(trim_data_filt(:,ckk),'-k');
        hold on;
    end
    %=====================================================================
    % write out the outputs
    scalars_out(1,:,ckk) = midpoints_red;
    scalars_out(2,:,ckk) = mitpsval_red;
    zlocs_out(:,ckk) = uint16(points_all);
    %orig_nsamples(ckk) = samples_arr;  
    
    %=====================================================================
    %uncompress the data to compare
    % for production just use the previous value rather than recomputing,
    % this is just code to put in a decompress function
    
    % need to have saved
    %     1) midpoints_red
    %     2) mitpsval_red
    %     3) scalepeak
    %     4) points_all
    
    %get the number of samples from the input
    %intrlength = size(trim_data,1);
    intrlength = intrlen;
    %samples_arr = single(1:intrlength)';
    %samples_arr = origcount;
    
    %get the value to scale to from the trace header
    %scalepeak = double(scalepeak);

    resmposenv = mposenv;
    %%resmposenv = double(interp1q(points_all,mitpsval_red,samples_arr));
    %resmposenv = double(makefastinterp1(double(points_all),double(mitpsval_red),double(samples_arr)));
    
    res_midpoints = midpoints;
    %%res_midpoints = double(interp1q(points_all,midpoints_red,samples_arr));
    %res_midpoints = double(makefastinterp1(double(points_all),double(midpoints_red),double(samples_arr)));
    
    scaleval_res = (res_midpoints - resmposenv)./scalepeak;
    scaleval_res(isnan(scaleval_res)) = 1;
    
    diffout = single((single(trim_data_filt(:,ckk)).*scaleval_res)+res_midpoints);
    
    if envlplot == 3
        figure(101); plot(trim_data(:,ckk),'-k');
        hold on;
        plot(diffout,'-r');
        figure(102); plot(abs((trim_data(:,ckk) - diffout)));
        %          hold on;
        %          figure(103); plot((abs((trim_data(:,ckk) - diffout)))./abs((trim_data(:,ckk))),'-g');
        
    end
    
    summoferror = sum(abs((trim_data(:,ckk) - diffout)));
    summofinput = sum(abs(trim_data(:,ckk)));
    if or(and( (summoferror/summofinput) > 0.0001, morecompress > 2), and( (summoferror/summofinput) > 0.01, morecompress == 2))
        bytes_to_store = ((size(mitpsval_red,1)*4) + ((size(points_all,1)*8) + 4));
        if morecompress < 3
            compression_ratio = (intrlength*4) / ((intrlength) + bytes_to_store);
        else
            compression_ratio = (intrlength*4) / ((intrlength*2) + bytes_to_store);
        end
        
        fprintf('trace %d ; %-10.8f percent error in reconstruction of dataset, one in %d, using %d bytes of scalars compression ratio %-5.2f \n',ckk,(summoferror/summofinput)*100,round((summofinput/summoferror)),bytes_to_store,compression_ratio);
    end
    %     if outplot == 1
    %         figure(50); plot(diffout);
    %         %hold on
    %         figure; plot(trim_data(:,ckk),'-g');
    %     end
end


end

