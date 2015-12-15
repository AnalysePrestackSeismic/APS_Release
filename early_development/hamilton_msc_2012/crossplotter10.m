function[] = crossplotter10(Yaxisfile,Xaxisfile, nt,ntwin,bytespersample, divisions, run_number, no_xlines,no_inlines, no_cross_bins,no_inline_bin)
fprintf('\n Reading in Data..\n');
%% READ IN DATA AND DETERMINE FILE SIZE/ NO. OF TRACES
fid1 = fopen(Yaxisfile);
fid2 = fopen(Xaxisfile);
fseek(fid1, 0, 'eof'); 
filesize1 = ftell(fid1) 
frewind(fid1); 
fseek(fid2, 0, 'eof');
filesize2 = ftell(fid2)
frewind(fid2);
ntrace1 = filesize1/(nt*bytespersample)
ntrace2 = filesize2/(nt*bytespersample)
fid4 = fopen(sprintf('%d_EEI_data',run_number),'a');
fid5 = fopen(sprintf('%d_flat_bk_trend_angle',run_number),'a');
%% Do the Crossplotting on binned data rather than whole slices.
% Number of time windows to look at
nloops = floor(nt/ntwin)
% Left over time windows
remainder = nt - (nloops*ntwin)
% Crossline Binsize
binsize = floor(no_xlines)/(no_cross_bins)
% Total Number of Samples per file
no_samples= (filesize1/4);
% Total Number of Samples per Y axis Slice
no_samples_slice = (filesize1/4)/nt;
% Total Number of Samples per X axis Slice
no_samples_slice1 = (filesize2/4)/nt;
% Counters
sum =  0;
hell = 1;
kickstarter= 0;
bin_jumper = no_samples_slice;
% hop jumps to the next inline
hop = (no_xlines);
binsize_in = 31;%no_inlines/no_inline_bin
starter = 0
count3 = 1
jump= 1;
sampleskip = floor(nt/divisions)-1;
count2 = 0;
inline_number_iterated = 0;
buffer2(:,1) = fread(fid1,no_samples_slice*466,'float32');   
buffer(:,1)= fread(fid2,no_samples_slice*466,'float32');
starter = 0;
finish = 0; 
nt_iterated = 0;
tic 
int_nan_finder = find(buffer(:,1) >1e29);
grad_nan_finder= find(buffer2(:,1) >1e29);                   
buffer(int_nan_finder) = NaN;
buffer(grad_nan_finder)= NaN;
buffer2(int_nan_finder) = NaN;
buffer2(grad_nan_finder)= NaN;
toc
numberofinlinebins = 0
time_counter=1

while finish ~= 1 
    % Reads in a Time window for whole slice of int and grad data
    jumper = 0;
    count2 = 0;  
    % Reads in all inlines on ntwin slices
    l = 1;
    % Starting location in terms of the number of inlines
    jumper = 0;
      for j = 1:ntwin                     
      for p = 1: binsize_in   
               % Location on current inline + location of current
               % xline + slice jumper
               a = starter+(p-1)*hop+((1)+(j-1)*bin_jumper);
               b = starter+((p-1)*hop+(((binsize)+(j-1)*bin_jumper)));     
               timeslice(l).graddata(jumper+(1:binsize),1) =  buffer2(a:b,1);
               timeslice(l).intdata(jumper+(1:binsize),1)  =  buffer(a:b,1);
               count2 = p;
               jumper = jumper + binsize; 
               
      end
      end
      if time_counter == 1
      tic 
      time_counter = 2
      end
    %numberofinlinebins = numberofinlinebins + 1
     
    if mod(b+1,no_xlines) == 0 
    starter = starter+(binsize)+1+(no_xlines);
    PCA(l).DEG
    inline_number_iterated = inline_number_iterated + 2
    if inline_number_iterated == (no_inlines-(binsize_in-1))
    starter = starter + (binsize_in*no_xlines)-(no_xlines)
    inline_number_iterated = 0;
    nt_iterated = nt_iterated+1
    toc
    tic
    end
    if nt_iterated == (466 - ntwin)
    finish = 1;
    end
    else
    starter = starter +2;
    end
                
    timeslice(l).intdata1(:,1)  = timeslice(l).intdata(isfinite(timeslice(l).intdata(:,1)),1);
    timeslice(l).graddata1(:,1) = timeslice(l).graddata(isfinite(timeslice(l).graddata(:,1)),1);         
    PCA(l).X = [timeslice(l).intdata1(:,1) timeslice(l).graddata1(:,1)];
    [PCA(l).COEFF,PCA(l).SCORE,PCA(l).latent]= princomp1(PCA(l).X);
    PCA(l).ANGLE = atan(-1*((PCA(l).COEFF(1,1)./PCA(l).COEFF(2,1))));
    EEI(1:length(timeslice(l).intdata)) = timeslice(l).intdata(:,1)*cos(PCA(l).ANGLE) + timeslice(l).graddata(:,1)*sin(PCA(l).ANGLE);
    PCA(l).DEG = 180*(PCA(l).ANGLE/pi);
    plotangle = PCA(l).DEG;
    princomp
        

clearvars timeslice
clearvars indexor
hell = 1;
fwrite(fid5,plotangle,'float32');
clearvars Array

end

end


