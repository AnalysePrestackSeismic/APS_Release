function[angle mean stddev] = crossplotter6(Yaxisfile,Xaxisfile, nt,ntwin,bytespersample, divisions, run_number, no_xlines,no_inlines, no_bins,no_inline_bin)

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


sampleskip = floor(nt/divisions)-1;

fid4 = fopen(sprintf('%d_EEI_data',run_number),'a');
fid5 = fopen(sprintf('%d_flat_bk_trend_angle',run_number),'a');

%% Video Writing Details

screen_size = get(0, 'ScreenSize');
writerObj=VideoWriter(sprintf('%d_PCA_crossplotmovie',run_number));
writerObj.FrameRate=24;
open(writerObj);

%% Do the Crossplotting on binned data rather than whole slices.

nloops = floor(nt/ntwin)
remainder = nt - floor(nloops*ntwin)

binsize = floor(no_xlines)/(no_bins)
no_samples= (filesize1/4);
no_samples_slice = (filesize1/4)/nt;
no_samples_slice1 = (filesize2/4)/nt;
sum =  0;
hell = 1;
kickstarter= 0;
x = -250:5:250;
bin_jumper = no_samples_slice;
count2 = 0;

hop = (no_xlines);

nloops2 = no_inlines/no_inline_bin
starter = 0
count3 = 1
x = -600:0.1:600

buffer(:,1) = fread(fid1,5*no_samples_slice*ntwin,'float32');
buffer2(:,1)= fread(fid2,5*no_samples_slice*ntwin,'float32');  
f(1:2*no_samples_slice) = 1e30;
fwrite(fid5,f,'float32');





for l = 1:nt - 4
    
    if l>1
    buffer(1:no_samples_slice) = [];
    buffer2(1:no_samples_slice) = [];
    buffer(4*no_samples_slice+1:5*no_samples_slice,1) = fread(fid1,no_samples_slice,'float32');
    buffer2(4*no_samples_slice+1:5*no_samples_slice,1,1)= fread(fid2,no_samples_slice,'float32');  
    end
    
    
        %Reads in all inlines on ntwin slices
               timeslice(l).graddata(:,1) =  buffer(:,1);
               timeslice(l).intdata(:,1)  =  buffer2(:,1);

        % Indexor 1 contains the full range of indices for the vector        
                     
int_nan_finder = find(buffer2(:,1) >1e29);
grad_nan_finder= find(buffer(:,1) >1e29);                   
buffer(int_nan_finder) = NaN;
buffer(grad_nan_finder)= NaN;
buffer2(int_nan_finder) = NaN;
buffer2(grad_nan_finder)= NaN;

%         for m = 1:length(timeslice(l).intdata)
%                         if  timeslice(l).intdata(m,1) > 1e29
%                                 timeslice(l).intdata(m,1) =  NaN;
%                                     timeslice(l).graddata(m,1) = NaN;
%             
%                         end
%                             if timeslice(l).graddata(m,1) > 1e29
%                                  timeslice(l).graddata(m,1) = NaN;
%                                     timeslice(l).intdata(m,1) =  NaN;
%                             end
%                       end
                   
          timeslice(l).intdata1(:,1)  = buffer2(isfinite(buffer2));
          timeslice(l).graddata1(:,1) = buffer(isfinite(buffer));
          
          
          
          
          hold on
          PCA(l).X = [timeslice(l).intdata1(:,1) timeslice(l).graddata1(:,1)];
          [PCA(l).COEFF,PCA(l).SCORE,PCA(l).latent]= princomp1(PCA(l).X);
          PCA(l).ANGLE = atan(-1*((PCA(l).COEFF(1,1)./PCA(l).COEFF(2,1))));
          EEI(1:length(timeslice(l).intdata)) = timeslice(l).intdata(:,1)*cos(PCA(l).ANGLE) + timeslice(l).graddata(:,1)*sin(PCA(l).ANGLE);
          PCA(l).DEG = 180*(PCA(l).ANGLE/pi);
          plotangle(1:length(timeslice(l).intdata)) = PCA(l).DEG;   

AngleArray(1:(no_samples_slice)) =PCA(l).DEG;
fwrite(fid5,AngleArray(:),'float32');

clearvars timeslice
clearvars indexor

end
fwrite(fid5,f,'float32');


% w = 1;
% 
% for g = 1:no_bins*no_inline_bin
% 
% for k = 1:nloops
% 
% bin_counter = no_bins*no_inline_bin   
% w = w + bin_counter;
% 
% plot(angle(g+(k-1)*bin_counter,1),k,'+')   
% set(gca,'YDir','reverse');
% hold on
% 
% 
% 
% end
% 
% end

end




