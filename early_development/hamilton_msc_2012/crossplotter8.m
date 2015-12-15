function[angle] = crossplotter8(Yaxisfile,Xaxisfile,Time_Volume, nt,ntwin,bytespersample, divisions, run_number, no_xlines,no_inlines, no_bins,no_inline_bin)

fprintf('\n Reading in Data..\n');
%% READ IN DATA AND DETERMINE FILE SIZE/ NO. OF TRACES

fid1 = fopen(Yaxisfile);
fid2 = fopen(Xaxisfile);
fid3 = fopen(Time_Volume);



fseek(fid1, 0, 'eof'); 
filesize1 = ftell(fid1) 
frewind(fid1); 

fseek(fid2, 0, 'eof');
filesize2 = ftell(fid2)
frewind(fid2);

fseek(fid3, 0, 'eof');
filesize3 = ftell(fid3)
frewind(fid3);

ntrace1 = filesize1/(nt*bytespersample)
ntrace2 = filesize2/(nt*bytespersample)
ntrace3 = filesize3/(nt*bytespersample)

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
remainder = nt - (nloops*ntwin)

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
jump= 1;


for k = 1:nloops
       
    no_samples_slice*ntwin
    buffer(:,1) = fread(fid1,no_samples_slice*ntwin,'float32');
    buffer2(:,1)= fread(fid2,no_samples_slice*ntwin,'float32');  
    count = 1;
    jumper = 0;
    count2 = 0;
    starter = 0;
    length(buffer);
    
    %Reads in all inlines on ntwin slices
    
    for l = 1:no_inline_bin   
       
        starter = starter + count2*no_xlines;    
    
        for i = 1:no_bins  
              jumper = 0;
         
            for j = 1:ntwin      
        
                
                for p = 1: nloops2
   
                    
                    
                a = starter+(p-1)*hop+(((i-1)*binsize+1)+(j-1)*bin_jumper);
                b = starter+((p-1)*hop+(((i*binsize)+(j-1)*bin_jumper)));
             
 
               timeslice(l).graddata(jumper+(1:binsize),1) =  buffer(a:b,1);
               timeslice(l).intdata(jumper+(1:binsize),1)  =  buffer2(a:b,1);
     
                count2 = p;
                jumper = jumper + binsize;
                
                % Need to change accordingly
                
                indices =   (starter+(p-1)*hop+(((count-1)*binsize+1)+(j-1)*bin_jumper):starter+((p-1)*hop+(((count*binsize)+(j-1)*bin_jumper))));           
                
                if j == 1 && p == 1
                    
                    indexor = indices;
                    
                else
                    
                    indexor = [ indexor indices];
                    
                end
               
                end
                
            end

        % Indexor 1 contains the full range of indices for the vector        
      
        
                    for m = 1:length(timeslice(l).intdata)
                        if  timeslice(l).intdata(m,1) > 1e29
                                timeslice(l).intdata(m,1) =  NaN;
                                    timeslice(l).graddata(m,1) = NaN;
            
                        end
                            if timeslice(l).graddata(m,1) > 1e29
                                 timeslice(l).graddata(m,1) = NaN;
                                    timeslice(l).intdata(m,1) =  NaN;
                            end
                    end                   
                    
          timeslice(l).intdata1(:,1)  = timeslice(l).intdata(isfinite(timeslice(l).intdata(:,1)),1);
          timeslice(l).graddata1(:,1) = timeslice(l).graddata(isfinite(timeslice(l).graddata(:,1)),1);
     
           
          PCA(l).X = [timeslice(l).intdata1(:,1) timeslice(l).graddata1(:,1)];
          [PCA(l).COEFF,PCA(l).SCORE,PCA(l).latent]= princomp1(PCA(l).X);
          PCA(l).ANGLE = atan(-1*((PCA(l).COEFF(1,1)./PCA(l).COEFF(2,1))));
          EEI(1:length(timeslice(l).intdata)) = timeslice(l).intdata(:,1)*cos(PCA(l).ANGLE) + timeslice(l).graddata(:,1)*sin(PCA(l).ANGLE);
          PCA(l).DEG = 180*(PCA(l).ANGLE/pi);
%         plotangle(1:length(timeslice(l).intdata)) = 0;
          angle_keeper(count3,1) = PCA(l).DEG;
          count3 = count3 + 1;
          count = count +1;
          plotangle = PCA(l).DEG;        
          
%                 subplot(4,1,1)
%                 plot(timeslice(l).intdata,timeslice(l).graddata,'+');
%                 hold on
%                 quiver(0,0,200*(PCA(l).COEFF(1,1)),200*(PCA(l).COEFF(2,1)),'r','linewidth',3);
%                 axis([-600 600 -600 600])   
%                 title(sprintf(' Inline %4.4f Bin %4.4f nloop %4.4f',l,i,k));
%                 hold off
%                
%                 subplot(4,1,2)
%                 n = hist(EEI,x);
%                 hist(EEI,x);
%                 title('EEI Histogram - Bin Size = 1')
%                 ph = get(gca,'children');
%                 vn = get(ph,'Vertices');
%                 vn(:,2) = vn(:,2) + 1;
%                 set(ph,'Vertices',vn);
%                 set(gca,'yscale','log')                      
%                 
%                 subplot(4,1,3)  
%                 m = hist(timeslice(l).intdata1,x);
%                 hist(timeslice(l).intdata1,x);
%                 ph = get(gca,'children');
%                 vn = get(ph,'Vertices');
%                 vn(:,2) = vn(:,2) + 1;
%                 set(ph,'Vertices',vn);
%                 set(gca,'yscale','log')       
%                 title('Intercept Histogram - Bin Size = 1')
%           
%                 subplot(4,1,4)   
%                 w = abs(m-n);       
%                 semilogy(x,w,'r+');
%                 axis([-250 250 0 10000])  
%                 title('Intercept Histogram - Bin Size = 1')       
%                 
%                 frame=getframe(figure(1));
%                 writeVideo(writerObj,frame);
%                 sum = length(timeslice(l).intdata1) + sum;
 

if hell == 1
%                       Array(:,1) = [indexor];
                        Array(1,1) = [plotangle];
%                       Array(:,3) = [EEI];
                        hell = hell +1;      

else
    
%                       B1 = Array(:,1);
%                       B2 = Array(:,1)    
%                       B3 = Array(:,3);    
                        Array(hell,1) = [plotangle];   
 
                        hell = hell +1;
end


clearvars timeslice
clearvars indexor

end
end
       
hell = 1;
% [ tmp ind] = sort(Array(:,1));Array = Array(ind,:);
fwrite(fid5,Array,'float32');
% fwrite(fid4,Array(:,3),'float32');
clearvars Array

end
% zerodata = zeros(round(remainder),1)
% fwrite(fid5,zerodata,'float32');
% fwrite(fid4,zerodata,'float32');
close(writerObj);
angle = angle_keeper;

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
% end
% 
% end

end




