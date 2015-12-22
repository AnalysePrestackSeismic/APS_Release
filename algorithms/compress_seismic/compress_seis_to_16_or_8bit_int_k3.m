function [ trim_data_filt, scalars_out, zlocs_out, tpgrad_orig, scalepeak,  orig_nsamples, percenterr ] = compress_seis_to_16_or_8bit_int_k3( trim_data, morecompress, envlplot )
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
% Defination: compress_seis_to_16_or_8bit_int compress the seismic into 16 bit
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

%trim_data = reshape(trim_data,[],floor(size(trim_data,2)/50));
no_of_traces = size(trim_data,2);
intrlen = size(trim_data,1);
intrlen_orig = intrlen;
orig_nsamples = uint16(intrlen);
percenterr = zeros(no_of_traces,1,'single');
% number of envelopes to try
nl = 20;
%outplot = 0;
firstpos = single(0);
endpos = single(0);
firstneg = single(0);
endneg = single(0);
%finaldroptri = 0.000001;
finaldroptri = 0.002;
%nanperwhtnoise  = 0.02;
diffpercentcut_orig = 0.2; %1.5; % this is the difference that it starts to drop at so bigger is more points dropped
diffpercentcut = diffpercentcut_orig;
diffpercentcut_update = 0.1;
diffpercentcut_max = 1.4;
nearfirst = 0;
% this should vary with trace length

if morecompress == 1
    tpgrad_orig = floor(intrlen/10);  % max number of scalars to store
    if tpgrad_orig < 55
        tpgrad_orig = 55;
    end
    extrapoints = floor(intrlen/200); % max number of extra pointas to hold to take care of points missed in derivative envelope
    if extrapoints < 8
        extrapoints = 8;
    end
    nanperwhtnoise  = 0.01;  % % to add to avoid cliping
    %perwhtnoise = 0.3;
    scalepeak = double(100);   % the max value to scale to
    smalldrop = 0.000000001;
elseif morecompress == 2
    tpgrad_orig = floor(intrlen/20);   % max number of scalars to store
    if tpgrad_orig < 30 % was 55
        tpgrad_orig = 30; % was 55
    end
    extrapoints = floor(intrlen/350); % max number of extra pointas to hold to take care of points missed in derivative envelope
    if extrapoints < 6
        extrapoints = 6;
    end
    extrapoints = 0;
    %fprintf('compression_space is %d bytes\n',((tpgrad_orig*10)+4));
    nanperwhtnoise  = 0.02; % 0.007;  % % to add to avoid cliping
    %perwhtnoise = 0.06;
    scalepeak = double(126); % the max value to scale to
    smalldrop = 0.000000000001;
    %diffpercentcut = floor(intrlen/600); %1.5; % this is the difference that it starts to drop at so bigger is more points dropped
    % this should vary with trace length
    % number of envelopes to try
    %nl = floor(intrlen/20);
else
    tpgrad_orig = floor(intrlen/50);   % max number of scalars to store
    if tpgrad_orig < 20
        tpgrad_orig = 20;
    end
    extrapoints = floor(intrlen/500); % max number of extra pointas to hold to take care of points missed in derivative envelope
    if extrapoints < 3
        extrapoints = 3;
    end
    nanperwhtnoise  = 0.01;  % % to add to avoid cliping
    %perwhtnoise = 0.01;
    scalepeak = double(28000); % the max value to scale to
    smalldrop = 0.000000001;
end

%to make the envelope tend to the value +/- 2000
%   Detailed explanation goes here
%%
midpoints_red_t = [ones(tpgrad_orig,1,'single') zeros(tpgrad_orig,1,'single')];
midpoints_red_torg = midpoints_red_t;
points_all_blk = ones(tpgrad_orig,1,'single');
points_all_blk_org = points_all_blk;

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
scaleval_res = ones(intrlen,1,'double');
res_midpoints = midblank;
scaleval = midblank;
%midpointslogic = logical(diffout);
midpoints = midblank;
itpslocs = single(1:intrlen)';
itpslocstest = itpslocs;
origcount = itpslocs;
origcountconst = origcount;
tmpupper = itpslocstest;
blankone = single(1);
blanklogic = true;
blankdbl = 0;
xgrad = blankone;
ygrad = blankone;
pointlogic = blanklogic;
top15grads = blankone;
points_all = blankone;
itpslocsin = blankone;
itpsval = blankone;
mid_lens = blankone;
midpoints_red  = blankone;
midpoints_red_grad  = blankone;
mitpslocsin  = blankone;
mitpsval = blankone;
mitpsval_red = blankone;
%
scalars_out = ones(2,tpgrad_orig,no_of_traces,'single');
zlocs_out = ones(tpgrad_orig,no_of_traces,'uint16');
tpgrad_origsm = tpgrad_orig - extrapoints;
%nanperwhtnoiseb = 0.0001;
addback = 1;
clipshift = 1.5;
trim_datacj = trim_data;
trim_data = trim_datacj*0;
%
%
%
for ckk = 1:no_of_traces
    %for ckk = 1227:1227
    %fprintf('trace %d \n',ckk);
    % find if the trace has zreros on the start or end and then remove
    % those as no need to compress them with this scheme
    if sum(trim_datacj(:,ckk)) ~= 0;
        tracezeros = (trim_datacj(:,ckk) ~= 0).*origcountconst;
        tracezerosout = tracezeros(tracezeros > 0);
        lastreal = tracezerosout(end);
        firstreal = tracezerosout(1);
        nooffzeros = firstreal - 1;
        nooflzeros = intrlen_orig - lastreal;
        noofrealsamps = (lastreal-firstreal)+1;
        trim_data(1:noofrealsamps,ckk) = trim_datacj(firstreal:lastreal,ckk);
        
        %set the other items to work for the reduced trace length
        intrlen = noofrealsamps;
        scaleval = midblank(1:noofrealsamps);
        midpoints = midblank(1:noofrealsamps);
        itpslocs = origcountconst(1:noofrealsamps);
        itpslocstest = origcountconst(1:noofrealsamps);
        origcount = origcountconst(1:noofrealsamps);
        diffpercentcut = diffpercentcut_orig;
        %tmpupper = origcount(1:noofrealsamps);
        
        % reset some variables
        midpoints_red_t = midpoints_red_torg;
        points_all_blk = points_all_blk_org;
        tpgrad = tpgrad_origsm;
        nearfirst = 0;
        lastdiff = noofrealsamps;
        % find the first and second derivatives of the max
        max_1st_deriv = diff(trim_data(1:noofrealsamps,ckk));
        max_2nd_deriv = diff(trim_data(1:noofrealsamps,ckk),2);
        
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
        % upper points
        itpsval = trim_data((diffsign < 0),ckk);
        itpslocsin = itpslocs(diffsign < 0);
        % lower points
        mitpsval = trim_data((mdiffsign > 0),ckk);
        mitpslocsin = itpslocs(mdiffsign > 0);
        
        %test to see if the envelope does not contain the start or end points
        % end of trace check
        if itpslocsin(end) > mitpslocsin(end)
            if itpslocsin(end) ~= intrlen
                itpslocsin = [itpslocsin; intrlen];
                %check to see if the real data extand outside of the last envelope
                %point
                if trim_data(lastreal,ckk) > itpsval(end)
                    itpsval = [itpsval; trim_data(lastreal,ckk)];
                else
                    itpsval = [itpsval; itpsval(end)];
                end
            end
            mitpslocsin = [mitpslocsin; intrlen];
            %check to see if the real data extand outside of the last envelope
            %point
            if trim_data(lastreal,ckk) < mitpsval(end)
                mitpsval = [mitpsval; trim_data(lastreal,ckk)];
            else
                mitpsval = [mitpsval; mitpsval(end)];
            end
        else
            if mitpslocsin(end) ~= intrlen
                mitpslocsin = [mitpslocsin; intrlen];
                %check to see if the real data extand outside of the last envelope
                %point
                if trim_data(lastreal,ckk) < mitpsval(end)
                    mitpsval = [mitpsval; trim_data(lastreal,ckk)];
                else
                    mitpsval = [mitpsval; mitpsval(end)];
                end
            end
            itpslocsin = [itpslocsin; intrlen];
            %check to see if the real data extand outside of the last envelope
            %point
            if trim_data(lastreal,ckk) > itpsval(end)
                itpsval = [itpsval; trim_data(lastreal,ckk)];
            else
                itpsval = [itpsval; itpsval(end)];
            end
        end
        
        %start of trace check
        if itpslocsin(1) < mitpslocsin(1)
            if itpslocsin(1) ~= 1
                itpslocsin = [1; itpslocsin];
                itpsval = [itpsval(1); itpsval];
            end
            mitpslocsin = [1; mitpslocsin];
            if trim_data(lastreal,ckk) < mitpsval(1)
                mitpsval = [trim_data(1,ckk); mitpsval];
            else
                mitpsval = [mitpsval(1); mitpsval];
            end
            if mitpsval(1) ==  itpsval(1)
                itpsval(1) = itpsval(2);
            end
        else
            if mitpslocsin(1) ~= 1
                mitpslocsin = [1; mitpslocsin];
                mitpsval = [mitpsval(1); mitpsval];
            end
            itpslocsin = [1; itpslocsin];
            if trim_data(lastreal,ckk) > itpsval(1)
                itpsval = [trim_data(1,ckk); itpsval];
            else
                itpsval = [itpsval(1); itpsval];
            end
            if mitpsval(1) ==  itpsval(1)
                mitpsval(1) = mitpsval(2);
            end
        end
        
        %store orginal points to compare at the end
        origuprloc = itpslocsin;
        origuprval = itpsval;
        origlowloc = mitpslocsin;
        origlowval = mitpsval;
        
        if envlplot == 4
            figure(44);
            plot(origlowloc,origlowval,'-or');
            hold on;
            plot(origuprloc,origuprval,'-og');
            plot(trim_data(1:noofrealsamps,ckk),'-k');
        end
        
        % make sure the start and end do not get lost
        firstpos = itpsval(1);
        endpos = itpsval(end);
        firstneg = mitpsval(1);
        endneg = mitpsval(end);
        
        % loop around the evelopes until dropped enough points ==========================================================================
        
        % test to remove points that are tiny and repeating, ie noise
        [itpslocsin, itpsval] = dropsmalldiff(itpslocsin,itpsval,finaldroptri);
        % test to remove points that are tiny and repeating, ie noise
        [mitpslocsin, mitpsval] = dropsmalldiff(mitpslocsin,mitpsval,finaldroptri);
        
        
        
        for envl = 1:nl
            %store orginal points to compare at the end
            origuprloc = itpslocsin;
            origuprval = itpsval;
            origlowloc = mitpslocsin;
            origlowval = mitpsval;
            %==========================================================================
            % upper envelope           
            % calc the envelope of the envelope for the upper envelope
            %[itpslocsin, itpsval] = upr_envlop(posenv, itpslocs,firstpos, endpos, intrlen, origcount, origuprloc, origuprval, itpslocstest);
            [itpslocsin, itpsval] = upr_envlop(itpsval, itpslocsin,firstpos, endpos, intrlen, origuprloc, origuprloc, origuprval, itpslocstest);
            
            % add small white noise to avoid division errors/mess
            itpsval = itpsval + abs(itpsval.*nanperwhtnoise);
            % final test to remove points that are tiny and repeating, ie noise
            %[itpslocsin, itpsval] = dropsmalldiff(itpslocsin,itpsval,finaldroptri);
            
            % interpolate to make the envelope only using fast linear interp
            posenv = double(makefastinterp1(double(itpslocsin),double(itpsval),double(origcount)));
            
            %==========================================================================
            % lower envelope           
            % calc the envelope of the envelope for the lower envelope
            
            [mitpslocsin, mitpsval] = lwr_envlop(mitpsval, mitpslocsin,firstneg, endneg, intrlen, origlowloc, origlowloc, origlowval, itpslocstest);
            %[mitpslocsin, mitpsval] = lwr_envlop(mposenv, itpslocs,firstneg, endneg, intrlen, origcount, origlowloc, origlowval, itpslocstest);
            
            %was adding some white noise
            mitpsval = mitpsval - abs(mitpsval.*nanperwhtnoise);
            % final test to remove points that are tiny and repeating, ie noise
            %[mitpslocsin, mitpsval] = dropsmalldiff(mitpslocsin,mitpsval,finaldroptri);
            
            % interpolate to make the envelope only using fast linear interp
            mposenv = double(makefastinterp1(double(mitpslocsin),double(mitpsval),double(origcount)));
            
            %=========================================================================
            if envlplot == -2
                figure(1105);
                plot(itpslocsin,itpsval,'-or');
                hold on;
                plot(mitpslocsin, mitpsval,'-or');
                %plot(updiflocs, updifval,'ob');
                plot(trim_data(1:noofrealsamps,ckk),'-k')
                plot(origuprloc,origuprval,'-g')
                plot(origlowloc,origlowval,'-g')
                hold off;
            end
            
            %=========================================================================
            % Take the difference between the 2nd and 1st envelope then find the envelope of that difference to find the key points to keep and add them back to the scalars
            % but try to make then centred round zero so the importance of
            % the misfit is maintained on upper and lower halves, so need
            % to remove or add the mid points to make a dc shift
            newnegvals =  mposenv(origlowloc)  - origlowval ; % difference to orginal lower values
            newposvals =  posenv(origuprloc) - origuprval ; % difference to the orginal upper values
            
            newnegvalsall = double(makefastinterp1(double(origlowloc),double(newnegvals),double(origcount)));
            newposvalsall = double(makefastinterp1(double(origuprloc),double(newposvals),double(origcount)));
            midpoints = newposvalsall+(((newposvalsall.*-1) - (newnegvalsall.*-1)).*0.5);
            newnegvalsall = newnegvalsall - midpoints;
            newposvalsall = newposvalsall - midpoints;
            newnegvals =  newnegvalsall(origlowloc);
            newposvals =  newposvalsall(origuprloc);
            % workout the distance between top and bottom envelopes and use it to threshold the envelope of the difference
            %totalenvdistup = posenv(origuprloc) - mposenv(origuprloc);
            %totalenvdistlow = posenv(origlowloc) - mposenv(origlowloc);
            %=========================================================================
            % Upper envelope calculate the envelope of the difference between the previous 2 envelopes the last variable 0 is upper 1 is lower envelope
            %[updiflocs, updifval] = diffenvlop(newposvals, origuprval,firstpos, endpos, intrlen, origuprloc, nanperwhtnoise, totalenvdistup, diffpercentcut, 0);
            [updiflocs, updifval] = diffenvlop(newposvals, origuprval,firstpos, endpos, intrlen, origuprloc, nanperwhtnoise, 0, diffpercentcut, 0);
            
            % add the two sets of points together
            [itpslocsin, itpsval] = addpoints(itpslocsin,itpsval,updiflocs,updifval,origcount);
itpslocsinplt = itpslocsin;
itpsvalplt = itpsval;            
            % add in the missed off points from the orginal envelope, but only those which are positive inflection
            [itpslocsin, itpsval] = outsidetest(origuprloc,origuprval,itpslocsin,itpsval,origcount,0,nanperwhtnoise);
            
            
            %=========================================================================
            % Lower envelope calculate the envelope of the difference between the previous 2 envelopes
            %[updiflocs, updifval] = diffenvlop(newnegvals, origlowval,firstneg, endneg, intrlen, origlowloc, nanperwhtnoise, totalenvdistlow, diffpercentcut, 1);
            [updiflocs, updifval] = diffenvlop(newnegvals, origlowval,firstneg, endneg, intrlen, origlowloc, nanperwhtnoise, 0, diffpercentcut, 1);
            
            % add the two sets of points together
            [mitpslocsin, mitpsval] = addpoints(mitpslocsin,mitpsval,updiflocs,updifval,origcount);
mitpslocsinplt = mitpslocsin;
mitpsvalplt = mitpsval;
            if envlplot == -1
                figure(1103);
                plot(itpslocsin,itpsval,'-or');
                hold on;
                plot(mitpslocsin, mitpsval,'or');
                %plot(updiflocs, updifval,'ob');
                plot(trim_data(1:noofrealsamps,ckk),'-k')
                %plot(updiflocs,updifval,'ob')
                plot(origuprloc,origuprval,'-g')
                plot(origlowloc,origlowval,'-g')
                hold off;
            end            
            
            % add in the missed off points from the orginal envelope, but only those which are positive inflection
            [mitpslocsin, mitpsval] = outsidetest(origlowloc,origlowval,mitpslocsin,mitpsval,origcount,1,nanperwhtnoise);
            
            if envlplot == -1
                figure(1102);
                plot(itpslocsin,itpsval,'-or');
                hold on;
                plot(mitpslocsin, mitpsval,'-or');
                %plot(updiflocs, updifval,'ob');
                plot(trim_data(1:noofrealsamps,ckk),'-k')
                plot(origuprloc,origuprval,'-g')
                plot(origlowloc,origlowval,'-g')
                hold off;
            end
            if envlplot == -3
                figure(2222);
                plot(origuprloc,origuprval,'-g')
                hold on;
                plot(origlowloc,origlowval,'-g')
                plot(mitpslocsinplt,mitpsvalplt,'-b')
                plot(itpslocsinplt,itpsvalplt,'-b')
                plot(itpslocsin,itpsval,'-or');
                plot(mitpslocsin, mitpsval,'-or');
                plot(updiflocs, updifval,'ob');
                plot(trim_data(1:noofrealsamps,ckk),'-k')
                plot(mposenv,'-m');
                plot(posenv,'-m');
                title(strcat(sprintf('loop no %d',envl),' green - orginal envelope; red new envelope; blue before fixing envelope'));
                hold off;
            end
            %=========================================================================
            
            
            % add in the missed off points from the orginal data, but only those which are positive inflection
            [mitpslocsin, mitpsval] = outsidetest_data(origcount,trim_data(1:noofrealsamps,ckk),mitpslocsin,mitpsval,origcount,1,(nanperwhtnoise*3));
            
            % add in the missed off points from the orginal data, but only those which are positive inflection
            [itpslocsin, itpsval] = outsidetest_data(origcount,trim_data(1:noofrealsamps,ckk),itpslocsin,itpsval,origcount,0,(nanperwhtnoise*3));
            
            %=========================================================================
            if envlplot == -2
                figure(1101);
                plot(itpslocsin,itpsval,'-or');
                hold on;
                plot(mitpslocsin, mitpsval,'-or');
                %plot(updiflocs, updifval,'ob');
                plot(trim_data(1:noofrealsamps,ckk),'-k')
                plot(origuprloc,origuprval,'-g')
                plot(origlowloc,origlowval,'-g')
                hold off;
            end
            
            %=========================================================================
            
            % decide if time to exit the loop
            if (size(mitpslocsin,1) + size(itpslocsin,1)) < floor(tpgrad);
                fprintf('under limit no of points %d ; no of points used %d\n',floor(tpgrad),(size(mitpslocsin,1) + size(itpslocsin,1)));
                break
            elseif  (size(mitpslocsin,1) + size(itpslocsin,1)) < floor(tpgrad)*1.3;
                if nearfirst == 1
                    diffpercentcut = diffpercentcut_orig*2;
                    nearfirst = 0;
                else
                    % if the difference between the last update and this one is small boost the update increase
                    if lastdiff - (size(mitpslocsin,1) + size(itpslocsin,1)) < (tpgrad*0.1)
                        diffpercentcut = diffpercentcut + (diffpercentcut_update*2);
                    else
                        diffpercentcut = diffpercentcut + diffpercentcut_update;
                    end
                end
                if diffpercentcut > diffpercentcut_max
                    diffpercentcut = diffpercentcut_max;
                end
                fprintf('near limit no of points; diffcut %f ; no of points used %d\n',diffpercentcut,(size(mitpslocsin,1) + size(itpslocsin,1)));
                
            else
                fprintf('limit no of points %d ;diffcut %f ; no of points used %d\n',floor(tpgrad),diffpercentcut,(size(mitpslocsin,1) + size(itpslocsin,1)));
                % if the difference between the last update and this one is small boost the update increase
                if lastdiff - (size(mitpslocsin,1) + size(itpslocsin,1)) < (tpgrad*0.25)
                    diffpercentcut = diffpercentcut + (diffpercentcut_update*2);
                else
                    diffpercentcut = diffpercentcut + diffpercentcut_update;
                end
                %diffpercentcut = diffpercentcut + diffpercentcut_update;
                if diffpercentcut > diffpercentcut_max
                    diffpercentcut = diffpercentcut_max;
                end
                nearfirst = 1;
            end
            lastdiff = size(mitpslocsin,1) + size(itpslocsin,1);
        end
        
        
        %==========================================================================
        if envlplot == -1
            figure(10001);
            plot(itpslocsin,itpsval,'-or');
            hold on;
            plot(mitpslocsin, mitpsval,'-or');
            %plot(updiflocs, updifval,'ob');
            plot(trim_data(1:noofrealsamps,ckk),'-k')
            %         plot(origuprloc,origuprval,'-g')
            %         plot(origlowloc,origlowval,'-g')
            hold off;
        end
        
        %reinterpolate all the points to be able to find the mid points
        mposenv = double(makefastinterp1(double(mitpslocsin),double(mitpsval),double(origcount)));
        posenv = double(makefastinterp1(double(itpslocsin),double(itpsval),double(origcount)));
        
        
        % calculate the mid points of the upper and lower envelopes from the
        % interpolated envelopes
        midpoints = posenv+(((posenv.*-1) - (mposenv.*-1)).*0.5);
        totaldist = (posenv - mposenv) + abs(midpoints);
        
        if envlplot == 10
            figure(5); plot(midpoints,'-b'); hold on
            %plot(midpoints_red_t(:,1),midpoints_red_t(:,2),'or');
            plot(posenv,'-r')
            plot(totaldist,'-m')
            %plot(midpoints_red_t(:,1),mposenv(midpoints_red_t(:,1)),'-or');
            plot(mposenv,'-g');
            plot(trim_data(1:noofrealsamps,ckk),'-k')
        end
        
        %decimate the mid points and then reinterpolate same as the
        %uncompression would do
        %if morecompress < 3
        % make the mid points for all the points on up and lower envelope
        itpslocs(itpslocsin) = itpsval;
        itpslocs(mitpslocsin) = mitpsval;
        points_all = itpslocstest(itpslocs ~= itpslocstest);
        
        tmpdiff = origcount;
        tmpdiff(itpslocsin) = [0;diff(itpsval)];
        tmpdiff(mitpslocsin) = [0;diff(mitpsval)];
        diff_all = tmpdiff(tmpdiff ~= itpslocstest);
        %else
        %    points_all = mitpslocsin;
        %end
        
        % make the sets of points to store this works out the size of
        % envelope plus distance from zero of the mid point
        midpoints_red = single(midpoints(points_all));
        mitpsval_red = single(mposenv(points_all));
        totaldist_red = single(totaldist(points_all)+abs(diff_all));
        
        
        %         figure(3000); plot(points_all,midpoints_red,'-ob'); hold on
        %         %plot(midpoints_red_t(:,1),midpoints_red_t(:,2),'or');
        %         plot(points_all,mitpsval_red,'-or')
        %         %plot(midpoints_red_t(:,1),mposenv(midpoints_red_t(:,1)),'-or');
        %         plot(points_all,totaldist_red,'-og');
        %         plot(trim_data(1:noofrealsamps,ckk),'-k')
        %         plot(posenv,'-r')
        %         plot(mposenv,'-g');
        
        
        
        
        
        midpoints = double(makefastinterp1(double(points_all),double(midpoints_red),double(origcount)));
        midpoints(isnan(midpoints)) = 0;
        mposenv = double(makefastinterp1(double(points_all),double(mitpsval_red),double(origcount)));
        mposenv(isnan(mposenv)) = 0;
        
        if envlplot == 3
            figure(97); plot(trim_data(1:noofrealsamps,ckk),'-k');
            hold on;
            plot(mposenv,'-r');
            plot((midpoints-(mposenv-midpoints)),'-g');
            %plot(posenv,'-m');
            plot(midpoints,'-b');
            %plot((midpoints - mposenv),'-m');
        end
        
        
        %Apply the scaling to the input data
        %scalepos = 2000 ./ maxenv;
        %scalepos = bsxfun(@rdivide,2000,posenv);
        scaleval = scalepeak./(midpoints - mposenv);
        scaleval(isnan(scaleval)) = 1;
        scaleval(isinf(scaleval)) = 1;
        
        if morecompress < 3
            trim_data_filt(firstreal:lastreal,ckk) = int8((trim_data(1:noofrealsamps,ckk)-midpoints).*scaleval);
        else
            trim_data_filt(firstreal:lastreal,ckk) = int16((trim_data(1:noofrealsamps,ckk)-midpoints).*scaleval);
        end
    else % traces are entirely zeros
        if morecompress < 3
            trim_data_filt(:,ckk) = int8(trim_data(:,ckk));
        else
            trim_data_filt(:,ckk) = int16(trim_data(:,ckk));
        end
        nooffzeros = 0;
        nooflzeros = 0;
        mposenv = midpoints;
    end
    if envlplot == 3
        figure(99); plot(trim_data_filt(:,ckk),'-k');
        figure(999); plot(trim_data(:,ckk),'-ok');
        hold on;
        plot(mposenv,'-r');
        plot((midpoints-(mposenv-midpoints)),'-g');
    end
    %=====================================================================
    % write out the outputs
    if size(midpoints_red,1) > tpgrad_orig
        scalars_out(1,1:tpgrad_orig,ckk) = midpoints_red(1:tpgrad_orig);
        fprintf('too many midpoints_red values at trace %d \n',ckk);
    else
        scalars_out(1,1:size(midpoints_red,1),ckk) = midpoints_red;
    end
    if size(mitpsval_red,1) > tpgrad_orig
        scalars_out(2,1:tpgrad_orig,ckk) = mitpsval_red(1:tpgrad_orig);
        fprintf('too many mitpsval_red values at trace %d \n',ckk);
    else
        scalars_out(2,1:size(mitpsval_red,1),ckk) = mitpsval_red;
    end
    % add on the first live sample numder to the points for decompression
    points_all = points_all+nooffzeros;
    if size(points_all,1) > tpgrad_orig
        zlocs_out(1:tpgrad_orig,ckk) = uint16(points_all(1:tpgrad_orig));
        fprintf('too many points_all values at trace %d \n',ckk);
    else
        zlocs_out(1:size(points_all,1),ckk) = uint16(points_all);
    end
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
    %intrlength = (intrlen_orig- tpgrad_orig);
    intrlength = intrlen_orig;
    %samples_arr = single(1:intrlength)';
    %samples_arr = origcount;
    
    %now remember to remove the scalars from the bottom of the trace
    %trim_data_filt(:,ckk) = trim_data_filt(1:intrlength,ckk);
    
    %get the value to scale to from the trace header
    %scalepeak = double(scalepeak);
    
    % read the scaling points and samples from the trace removing dupliates
    %     points_all = single(zlocs_out(:,ckk));
    %     midpoints_red = single(scalars_out(1,:,ckk));
    %     midpoints_red = midpoints_red([1;diff(points_all)] > 0);
    %     mitpsval_red = single(scalars_out(2,:,ckk));
    %     mitpsval_red = mitpsval_red([1;diff(points_all)] > 0);
    %     points_all = points_all([1;diff(points_all)] > 0);
    
    % before interpolating need to add extra points at the start and end of
    % the trace as the scalars can be anywhere along the trace
    %     if(points_all(1) ~= 1)
    %         points_all = [1;points_all];
    %         mitpsval_red = [mitpsval_red(1);mitpsval_red];
    %         midpoints_red = [midpoints_red(1);midpoints_red];
    %     end
    
    resmposenv = mposenv;
    %%resmposenv = double(interp1q(points_all,mitpsval_red,samples_arr));
    %resmposenv = double(makefastinterp1(double(points_all),double(mitpsval_red),double(samples_arr)));
    
    res_midpoints = midpoints;
    %%res_midpoints = double(interp1q(points_all,midpoints_red,samples_arr));
    %res_midpoints = double(makefastinterp1(double(points_all),double(midpoints_red),double(samples_arr)));
    
    
    %scaleval_restmp = (res_midpoints - resmposenv)./scalepeak;
    scaleval_restmp = (midpoints - mposenv)./scalepeak;
    scaleval_restmp(isnan(scaleval_res)) = 1;
    % this expansion is just for this code, the decompression needs to just add the points into a set array of blanks or ones
    scaleval_res = [ones(nooffzeros,1,'double');scaleval_restmp;ones(nooflzeros,1,'double')];
    res_midpoints = [zeros(nooffzeros,1,'double');midpoints;zeros(nooflzeros,1,'double')];
    
    diffout = single((single(trim_data_filt(:,ckk)).*scaleval_res)+res_midpoints);
    
    if envlplot > 2
        figure(101); plot([zeros(nooffzeros,1,'single');trim_data(1:noofrealsamps,ckk);zeros(nooflzeros,1,'single')],'-k');
        hold on;
        plot(diffout,'-r');
        hold off;
        %figure(103); hist(single(trim_data_filt(:,ckk)),255);
        %figure(103); hist(single(trim_data_filt(:,ckk)),65536);
        %       figure(102); plot((abs(double([zeros(nooffzeros,1,'single');trim_data(1:noofrealsamps,ckk);zeros(nooflzeros,1,'single')]) - double(diffout)))./abs(double([zeros(nooffzeros,1,'single');trim_data(1:noofrealsamps,ckk);zeros(nooflzeros,1,'single')])),'-k');
        %          hold on;
        %          figure(103); plot((abs((trim_data(:,ckk) - diffout)))./abs((trim_data(:,ckk))),'-g');
        
    end
    if sum(trim_data(:,ckk)) ~= 0;
        % added a 2% skip on error check on start and end of trace incase it is
        % a truncated trace which needs error checking with input samples
        % (which will take time) and they do not matter vey much as just the
        % ends of the traces
        errsskip = ceil(noofrealsamps*0.02);
        percenterr(ckk) = single((abs(sum(abs(double(trim_data((1+errsskip):(noofrealsamps-errsskip),ckk)))) - sum( abs( double(diffout((firstreal+errsskip):(lastreal-errsskip))))))/sum(abs(double(trim_data((1+errsskip):(noofrealsamps-errsskip),ckk)))))*100);
        
        %percenterr(ckk) = single((abs(sum(abs(double(trim_data(1:noofrealsamps,ckk)))) - sum( abs( double(diffout))))/sum(abs(double(trim_data(1:noofrealsamps,ckk)))))*100);
        
        %     summoferror = sum(abs((trim_data(:,ckk) - diffout)));
        %     summofinput = sum(abs(trim_data(:,ckk)));
        %%if or(and( (summoferror/summofinput) > 0.0001, morecompress > 2), and( (summoferror/summofinput) > 0.01, morecompress == 2))
        if or(and( percenterr(ckk) > 0.25, morecompress > 2), and( percenterr(ckk) > 0.0000005, morecompress <= 2))
        %if or(and( percenterr(ckk) > 0.25, morecompress > 2), and( percenterr(ckk) > 0.005, morecompress <= 2))
            bytes_to_store = ((size(scalars_out,2)*8) + ((size(zlocs_out,1)*2) + 4));
            if morecompress < 3
                compression_ratio = (intrlength*4) / ((intrlength) + bytes_to_store);
            else
                compression_ratio = (intrlength*4) / ((intrlength*2) + bytes_to_store);
            end
            
            fprintf('trace %d ; %-10.8f percent error in reconstruction of dataset, one in %d, using %d bytes of scalars compression ratio %-5.2f \n',ckk,percenterr(ckk),round((100/percenterr(ckk))),bytes_to_store,compression_ratio);
        end
    end
    %
end

%==========================================================================
    function [fuprlocs, fuprsval] = upr_envlop(upr_in_arr, arr_in_locs,ffirstpos, fendpos, trc_len, forigcount, f_origuprloc, f_origuprval, fitpslocstest)
        
        % calc the envelope of the envelope
        fmax_1st_deriv = diff(upr_in_arr);
        fmax_2nd_deriv = diff(upr_in_arr,2);
        fsign_1_deriv = sign(fmax_1st_deriv);
        fsign_1_deriv(fsign_1_deriv == 0) = 1;
        fdiffsign = diff(fsign_1_deriv);
        fdiffsign(sign(fmax_2nd_deriv) > 0) = 0;
        fdiffsign = [1;fdiffsign];
        
        fuprsval = upr_in_arr(fdiffsign < 0);
        fuprlocs = arr_in_locs(fdiffsign < 0);
        % make sure the start and end do not get lost
        if fuprlocs(end) ~= trc_len
            fuprlocs = [fuprlocs; trc_len];
            fuprsval = [fuprsval; fendpos];
        end
        if fuprlocs(1) ~= 1
            fuprlocs = [1; fuprlocs];
            fuprsval = [ffirstpos; fuprsval];
        end
        
        %         % interpolate to make the envelope only using fast linear interp
        %         %posenv = double(interp1q(itpslocsin,itpsval,origcount));
        %         upr_in_arr = double(makefastinterp1(double(fuprlocs),double(fuprsval),double(forigcount)));
        %         %==========================================================================
        %         %test to see if any points have been chopped off in the envelope
        %         if (size(forigcount,1) == size(f_origuprloc,1))
        %             fuprsval = f_origuprval(upr_in_arr <= f_origuprval);
        %             fuprlocs = forigcount(upr_in_arr <= f_origuprval);
        %         else
        %             fnewposvals =  upr_in_arr(f_origuprloc);
        %             tmpsampseq = forigcount;
        %             tmpsampseq(fuprlocs) = fuprsval;
        %             tmpsampseq(f_origuprloc(fnewposvals < f_origuprval)) = f_origuprval(fnewposvals < f_origuprval);
        %             fuprlocs = fitpslocstest(tmpsampseq ~= fitpslocstest);
        %             fuprsval = tmpsampseq(tmpsampseq ~= fitpslocstest);
        %         end
    end

    function [flowlocs, flowval] = lwr_envlop(low_in_arr, l_arr_in_locs,l_firstneg, l_endneg, l_trc_len, florigcount, f_origlowloc, f_origlowval, flitpslocstest)
        % calc the envelope of the envelope
        flmax_1st_deriv = diff(low_in_arr);
        flmax_2nd_deriv = diff(low_in_arr,2);
        flsign_1_deriv = sign(flmax_1st_deriv);
        flsign_1_deriv(flsign_1_deriv == 0) = 1;
        fmdiffsign = diff(flsign_1_deriv);
        fmdiffsign(sign(flmax_2nd_deriv) <= 0) = 0;
        fmdiffsign = [1;fmdiffsign];
        flowval = low_in_arr(fmdiffsign > 0);
        flowlocs = l_arr_in_locs(fmdiffsign > 0);
        
        % make sure the start and end do not get lost
        if flowlocs(end) ~= l_trc_len
            flowlocs = [flowlocs; l_trc_len];
            flowval = [flowval; l_endneg];
        end
        if flowlocs(1) ~= 1
            flowlocs = [1; flowlocs];
            flowval = [l_firstneg; flowval];
        end
        
        %         % interpolate to make the envelope only using fast linear interp
        %         %mposenv = double(interp1q(mitpslocsin,mitpsval,origcount));
        %         low_in_arr = double(makefastinterp1(double(flowlocs),double(flowval),double(florigcount)));
        %         %==========================================================================
        %         %test to see if any points have been chopped off in the envelope
        %         if (size(florigcount,1) == size(f_origlowloc,1))
        %             flowval = f_origlowval(low_in_arr >= f_origlowval);
        %             flowlocs = florigcount(low_in_arr >= f_origlowval);
        %         else
        %             fnewlowvals =  low_in_arr(f_origlowloc);
        %             tmpsamplseq = florigcount;
        %             tmpsamplseq(flowlocs) = flowval;
        %             tmpsamplseq(f_origlowloc(fnewlowvals > f_origlowval)) = f_origlowval(fnewlowvals > f_origlowval);
        %             flowlocs = flitpslocstest(tmpsamplseq ~= flitpslocstest);
        %             flowval = tmpsamplseq(tmpsamplseq ~= flitpslocstest);
        %         end
    end

    function [fupdiflocs, fupdifval] = diffenvlop(fnewposvals, foriguprval,d_firstdif, d_enddif, d_trc_len, foriguprloc, fnanperwhtnoiseb, ~, fdiffpercentcut, edmode)
        % calc the envelope of the difference
        dmax_1st_deriv = diff(fnewposvals);
        dmax_2nd_deriv = diff(fnewposvals,2);
        dsign_1_deriv = sign(dmax_1st_deriv);
        dsign_1_deriv(dsign_1_deriv == 0) = 1;
        ddiffsign = diff(dsign_1_deriv);
        
        seqno = 1:1:size(foriguprloc,1);
        uplocs = seqno([-1;ddiffsign;-1] < 0);
        downlocs = seqno([1;ddiffsign;1] > 0);
        if size(downlocs,2) > size(uplocs,2)
            downlocs = downlocs(:,1:size(uplocs,2));
        elseif size(uplocs,2) > size(downlocs,2)
            uplocs = uplocs(:,1:size(downlocs,2));
        end
        
        if edmode == 0 % upper envelope
            %figure(3000);
            %plot(foriguprloc, foriguprval,'-or');
            %hold on;
            %plot(foriguprloc, fnewposvals,'-ok');
            
            
            fupdifval = foriguprval(uplocs(  (fnewposvals(uplocs).*[1;diff(foriguprloc(downlocs))])  > (fdiffpercentcut*mean(fnewposvals(uplocs).*[1;diff(foriguprloc(downlocs))])) )) ;

            %fupdiflocs = foriguprloc(uplocs(  (fnewposvals(uplocs).*[1;diff(foriguprloc(downlocs))])  > (fdiffpercentcut*mean(fnewposvals(uplocs).*[1;diff(foriguprloc(downlocs))])) )) ;            
            %plot(fupdiflocs, fupdifval,'ob');
            
            % to do still fupdifval = fupdifval + abs(fnewposvals(([-1;ddiffsign;-1] < 0 &  fnewposvals  > fdiffpercentcut*mean(fnewposvals([-1;ddiffsign;-1] < 0)) )).*fnanperwhtnoiseb);
            %fupdifval = fupdifval + abs(fnewposvals(uplocs(  (fnewposvals(uplocs).*[1;diff(foriguprloc(downlocs))])  > (fdiffpercentcut*mean(fnewposvals(uplocs).*[1;diff(foriguprloc(downlocs))])) ))    ).*fnanperwhtnoiseb;
            fupdifval = fupdifval + abs(fupdifval).*fnanperwhtnoiseb;
            fupdiflocs = foriguprloc(uplocs(  (fnewposvals(uplocs).*[1;diff(foriguprloc(downlocs))])  > (fdiffpercentcut*mean(fnewposvals(uplocs).*[1;diff(foriguprloc(downlocs))])) )) ;

            %plot(fupdiflocs, fupdifval,'og');
            %hold off;
        else % lower envelope
            %figure(3001);
            %plot(foriguprloc, foriguprval,'-or');
            %hold on;
            %plot(foriguprloc, fnewposvals,'-ok');
            
            fupdifval = foriguprval(downlocs(  (fnewposvals(downlocs).*[1;diff(foriguprloc(uplocs))])  < (fdiffpercentcut*mean(fnewposvals(downlocs).*[1;diff(foriguprloc(uplocs))])) )) ;
            %fupdiflocs = foriguprloc(downlocs(  (fnewposvals(downlocs).*[1;diff(foriguprloc(uplocs))])  < (fdiffpercentcut*mean(fnewposvals(downlocs).*[1;diff(foriguprloc(uplocs))])) )) ;
            %plot(fupdiflocs, fupdifval,'ob');
            
            % old to do still fupdifval = fupdifval - abs(fnewposvals(([-1;ddiffsign;-1] > 0 &  fnewposvals  < fdiffpercentcut*mean(fnewposvals([-1;ddiffsign;-1] > 0)) )).*fnanperwhtnoiseb);
            % this adds part of the differrence back
            %fupdifval = fupdifval - abs(fnewposvals(downlocs(  (fnewposvals(downlocs).*[1;diff(foriguprloc(uplocs))])  < (fdiffpercentcut*mean(fnewposvals(downlocs).*[1;diff(foriguprloc(uplocs))])) )) ).*fnanperwhtnoiseb;
            % this adds part of the orginal value back
            fupdifval = fupdifval - abs(fupdifval).*fnanperwhtnoiseb;
            fupdiflocs = foriguprloc(downlocs(  (fnewposvals(downlocs).*[1;diff(foriguprloc(uplocs))])  < (fdiffpercentcut*mean(fnewposvals(downlocs).*[1;diff(foriguprloc(uplocs))])) )) ;
            
            %plot(fupdiflocs, fupdifval,'og');
            %hold off;
        end
        % make sure the start and end do not get lost
        if isempty(fupdiflocs) ~= 1
            if fupdiflocs(end) ~= d_trc_len
                fupdiflocs = [fupdiflocs; d_trc_len];
                fupdifval = [fupdifval; d_enddif];
            end
            if fupdiflocs(1) ~= 1
                fupdiflocs = [1; fupdiflocs];
                fupdifval = [d_firstdif; fupdifval];
            end
        end
    end

    function [origlocs, origvals] = addpoints(origlocs,origvals,newlocs,newvals,fullorigcount)
        % add the two sets of points together
        if isempty(newlocs) ~= 1
            tmpmerge = fullorigcount;
            tmpmerge(origlocs) = origvals;  % populate times sequence with the orginal point values
            tmpmerge(newlocs) = newvals;  % populate times sequence with the new point values
            origlocs = fullorigcount(tmpmerge ~= fullorigcount);
            origvals = tmpmerge(tmpmerge ~= fullorigcount);
        end
    end

    function [toriglocs, torigvals] = outsidetest(toriglocs,torigvals,tnewlocs,tnewvals,florigcountb,ulmode,whitemul)
        % toriglocs,torigvals - these are the locactions and values of the previous envelope
        % tnewlocs,tnewvals - these are the locactions and values of the current envelope
        %test to see if any points have been chopped off in the envelope
        %florigcountb - this is the points of the whoel trace or just the sequenctial array which is the size of the orginal envelope
        newinterp = double(makefastinterp1(double(tnewlocs),double(tnewvals),double(florigcountb)));
        if ulmode == 0 % upper envelope
            %         if (size(florigcountb,1) == size(toriglocs,1))
            %             % need to finsih this bit
            %             newinterp(([-1;diff(torigvals,2);-1]<0) & (newinterp < torigvals)) = torigvals(([-1;diff(torigvals,2);-1]<0) & (newinterp < torigvals));
            %             plot(toriglocs,torigvals,'-g');
            %             %make logical index into the interpolated array
            %             toriglocs(tnewlocs) = 99999;
            %             toriglocs = toriglocs((toriglocs==99999)| ([-1;diff(torigvals,2);-1]<0) & (newinterp > torigvals));
            %             torigvals = newinterp((florigcountb == tnewlocs  )| ([-1;diff(torigvals,2);-1]<0) & (newinterp > torigvals));
            %         else
            dtmpupper = florigcountb;
            dtmpupper(tnewlocs) = tnewvals;
            %            dtmpupper(toriglocs(([-1;diff(torigvals,2);-1]<0) & (newinterp(toriglocs) < torigvals))) = torigvals(([-1;diff(torigvals,2);-1]<0) & (newinterp(toriglocs) < torigvals))+abs((torigvals(([-1;diff(torigvals,2);-1]<0) & (newinterp(toriglocs) < torigvals)).*whitemul));
            dtmpupper(toriglocs(( [diff(([-1;diff(torigvals,1)]./[-1;diff(toriglocs,1)]),1);-1] < 0) & (newinterp(toriglocs) < torigvals))) = torigvals(( [diff(([-1;diff(torigvals,1)]./[-1;diff(toriglocs,1)]),1);-1] < 0) & (newinterp(toriglocs) < torigvals))+abs((torigvals(( [diff(([-1;diff(torigvals,1)]./[-1;diff(toriglocs,1)]),1);-1] < 0) & (newinterp(toriglocs) < torigvals)).*whitemul));
            
            toriglocs = florigcountb(dtmpupper ~= florigcountb);
            torigvals = dtmpupper(dtmpupper ~= florigcountb);
            %         end
        else % lower envelope
            dtmpupper = florigcountb;
            dtmpupper(tnewlocs) = tnewvals;
            %             figure(3001);
            %             plot(toriglocs, torigvals,'-or');
            %             hold on;
            %             plot(toriglocs, newinterp(toriglocs),'ob');
            %             plot(toriglocs([-1;diff(torigvals,2);-1]>0),torigvals([-1;diff(torigvals,2);-1]>0),'og');
            %             hold off;
            %dtmpupper(toriglocs(([-1;diff(torigvals,2);-1]>0) & (newinterp(toriglocs) > torigvals))) = torigvals(([-1;diff(torigvals,2);-1]>0) & (newinterp(toriglocs) > torigvals))-abs((torigvals(([-1;diff(torigvals,2);-1]>0) & (newinterp(toriglocs) > torigvals)).*whitemul));
            dtmpupper(toriglocs(( [diff(([-1;diff(torigvals,1)]./[-1;diff(toriglocs,1)]),1);-1] > 0) & (newinterp(toriglocs) > torigvals))) = torigvals(( [diff(([-1;diff(torigvals,1)]./[-1;diff(toriglocs,1)]),1);-1] > 0) & (newinterp(toriglocs) > torigvals))-abs((torigvals(( [diff(([-1;diff(torigvals,1)]./[-1;diff(toriglocs,1)]),1);-1] > 0) & (newinterp(toriglocs) > torigvals)).*whitemul));
            toriglocs = florigcountb(dtmpupper ~= florigcountb);
            torigvals = dtmpupper(dtmpupper ~= florigcountb);
        end
    end

    function [toriglocs, torigvals] = dropsmalldiff(inlocs,invals,dropthresh)
        % final test to remove points that are tiny and repeating, ie noise
        %should really be a trim mean to drop the extreme points, but too
        %slow and the limit is set really really small
        %keeplocstmp = (inlocs -[0;inlocs(1:(end-1))]).*abs(invals -[0;invals(1:(end-1))]) + (abs(([inlocs(2:end);0]-inlocs)).*abs(invals -[invals(2:(end));0]));
        keeplocstmp = (inlocs -[0;inlocs(1:(end-1))]).*abs(invals -[0;invals(1:(end-1))]) + (abs(([inlocs(2:end);0]-inlocs)).*abs(invals -[invals(2:(end));0]));
        keeplocs =  keeplocstmp > (mean(keeplocstmp) * dropthresh);
        keeplocs(1:2) = 1;
        keeplocs(end) = 1;
        toriglocs = inlocs(keeplocs);
        torigvals = invals(keeplocs);
    end

    function [toriglocs, torigvals] = outsidetest_data(toriglocs,torigvals,tnewlocs,tnewvals,florigcountb,ulmode,whitemul)
        % toriglocs,torigvals - these are the locactions and values of the previous envelope
        % tnewlocs,tnewvals - these are the locactions and values of the current envelope
        %test to see if any points have been chopped off in the envelope
        %florigcountb - this is the points of the whoel trace or just the sequenctial array which is the size of the orginal envelope
        newinterp = double(makefastinterp1(double(tnewlocs),double(tnewvals),double(florigcountb)));
        if ulmode == 0 % upper envelope
            dtmpupper = florigcountb;
            dtmpupper(tnewlocs) = tnewvals;
            %            dtmpupper(toriglocs(([-1;diff(torigvals,2);-1]<0) & (newinterp(toriglocs) < torigvals))) = torigvals(([-1;diff(torigvals,2);-1]<0) & (newinterp(toriglocs) < torigvals))+abs((torigvals(([-1;diff(torigvals,2);-1]<0) & (newinterp(toriglocs) < torigvals)).*whitemul));
            dtmpupper(toriglocs(( [diff(([-1;diff(torigvals,1)]./[-1;diff(toriglocs,1)]),1);-1] < 0) & (newinterp(toriglocs) < torigvals))) = torigvals(( [diff(([-1;diff(torigvals,1)]./[-1;diff(toriglocs,1)]),1);-1] < 0) & (newinterp(toriglocs) < torigvals))+abs((torigvals(( [diff(([-1;diff(torigvals,1)]./[-1;diff(toriglocs,1)]),1);-1] < 0) & (newinterp(toriglocs) < torigvals)).*whitemul));
            
            toriglocs = florigcountb(dtmpupper ~= florigcountb);
            torigvals = dtmpupper(dtmpupper ~= florigcountb);
            %         end
        else % lower envelope
            dtmpupper = florigcountb;
            dtmpupper(tnewlocs) = tnewvals;
            %             figure(3001);
            %             plot(toriglocs, torigvals,'-or');
            %             hold on;
            %             plot(toriglocs, newinterp(toriglocs),'ob');
            dtmpupper(toriglocs(( [diff(([-1;diff(torigvals,1)]./[-1;diff(toriglocs,1)]),1);-1] > 0) & (newinterp(toriglocs) > torigvals))) = torigvals(( [diff(([-1;diff(torigvals,1)]./[-1;diff(toriglocs,1)]),1);-1] > 0) & (newinterp(toriglocs) > torigvals))-abs((torigvals(( [diff(([-1;diff(torigvals,1)]./[-1;diff(toriglocs,1)]),1);-1] > 0) & (newinterp(toriglocs) > torigvals)).*whitemul));
            toriglocs = florigcountb(dtmpupper ~= florigcountb);
            torigvals = dtmpupper(dtmpupper ~= florigcountb);
        end
    end
end

