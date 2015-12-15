function[] = crossplotter5(Yaxisfile,Xaxisfile,Time_Volume, nt,ntwin,bytespersample, divisions, run_number, no_xlines,no_inlines, no_bins)

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
remainder = no_xlines-(no_bins*binsize);
no_samples= (filesize1/4);
no_samples_slice = (filesize1/4)/nt;
no_samples_slice1 = (filesize2/4)/nt;
sum =  0;
hell = 1;
kickstarter= 0;
x = -600:10:600;

for k = 1:nloops
    
    
    buffer(:,1) = fread(fid1,no_samples_slice*ntwin,'float32');
    buffer2(:,1)= fread(fid2,no_samples_slice*ntwin,'float32');
    
    length(buffer)
    
    count = 1;

    %Reads in all inlines on ntwin slices
    for l = 1:no_inlines     
        %Reads in one bin on all slices on 1 inline

        for i = 1:no_bins
            %Reads in 1 Bin on all slices in window of 1 inline         

 
            for j = 1:ntwin       
               
               
                jumper = (j-1)*binsize;
                bin_jumper = no_samples_slice;
                                        
                timeslice(i).graddata(jumper+(1:binsize),1) =  buffer(kickstarter+(((count-1)*binsize+1)+(j-1)*bin_jumper):kickstarter+(((count*binsize)+(j-1)*bin_jumper)),1);
                timeslice(i).intdata(jumper+(1:binsize),1)  =  buffer2(kickstarter+(((count-1)*binsize+1)+(j-1)*bin_jumper):kickstarter+(((count*binsize)+(j-1)*bin_jumper)),1);
                
                indices =   (kickstarter+((count-1)*binsize+1)+(j-1)*bin_jumper):kickstarter+(((count*binsize)+(j-1)*bin_jumper));           
                
                if j == 1
                    
                    indexor = indices;
                    
                else
                    
                    indexor = [ indexor indices];
                    
                end
               
           end                             
           % Indexor 1 contains the full range of indices for the vector
           % plotangle              
            count = count + 1;
                    for m = 1:binsize*ntwin
                        if  timeslice(i).intdata(m,1) > 1e29
                                timeslice(i).intdata(m,i) =  NaN;
                                    timeslice(i).graddata(m,i) = NaN;
            
                        end
                            if timeslice(i).graddata(m,1) > 1e29
                                 timeslice(i).graddata(m,1) = NaN;
                                    timeslice(i).intdata(m,1) =  NaN;
                            end
                    end
                    
            
      
       
          

           
           
           timeslice(i).intdata1(:,1)  = timeslice(i).intdata(isfinite(timeslice(i).intdata(:,1)),1);
           timeslice(i).graddata1(:,1) = timeslice(i).graddata(isfinite(timeslice(i).graddata(:,1)),1);
            
            
            
            
            PCA(i).X = [timeslice(i).intdata1(:,1) timeslice(i).graddata1(:,1)];
            [PCA(i).COEFF,PCA(i).SCORE,PCA(i).latent]= princomp1(PCA(i).X);
            PCA(i).ANGLE = atan(-1*((PCA(i).COEFF(1,1)./PCA(i).COEFF(2,1))));
            EEI(1:length(timeslice(i).intdata)) = timeslice(i).intdata(:,1)*cos(PCA(i).ANGLE) + timeslice(i).graddata(:,1)*sin(PCA(i).ANGLE);
            PCA(i).DEG = 180*(PCA(i).ANGLE/pi);
            plotangle(1:length(timeslice(i).intdata)) = PCA(i).DEG;
            
                subplot(4,1,1)
                
                plot(timeslice(i).intdata,timeslice(i).graddata,'+');
                hold on
                quiver(0,0,200*(PCA(i).COEFF(1,1)),200*(PCA(i).COEFF(2,1)),'r','linewidth',3);
                axis([-600 600 -600 600])   
                title(sprintf(' Inline %4.4f Bin %4.4f nloop %4.4f',l,i,k));
                hold off
               
                subplot(4,1,2)
                n = hist(EEI,x);
                hist(EEI,x);
                title('EEI Histogram - Bin Size = 1')
                ph = get(gca,'children');
                vn = get(ph,'Vertices');
                vn(:,2) = vn(:,2) + 1;
                set(ph,'Vertices',vn);
                set(gca,'yscale','log')       
                               
                
                
                
                subplot(4,1,3)  
                m = hist(timeslice(i).intdata1,x);
                hist(timeslice(i).intdata1,x);
                ph = get(gca,'children');
                vn = get(ph,'Vertices');
                vn(:,2) = vn(:,2) + 1;
                set(ph,'Vertices',vn);
                set(gca,'yscale','log')       
                title('Intercept Histogram - Bin Size = 1')
          
                            
                
                subplot(4,1,4)   
                w = abs(m-n);       
                semilogy(x,w,'r+');
                axis([-600 600 0 1000])  
                
                
                
                
%                 ph = get(gca,'children');
%                 vn = get(ph,'Vertices');
%                 vn(:,2) = vn(:,2) + 1;
%                 set(ph,'Vertices',vn);
%                 set(gca,'yscale','log');      
                title('Intercept Histogram - Bin Size = 1');
             
                
                
                
                frame=getframe(figure(1));
                writeVideo(writerObj,frame);
                sum = length(timeslice(i).intdata1) + sum;
                       
            if hell == 1
               
                       Array(:,1) = [indexor];
                       Array(:,2) = [plotangle];
                       Array(:,3) = [EEI];
            else
               
                     B1 = Array(:,1);
                     B2 = Array(:,2);    
                     B3 = Array(:,3);    
                     Array = [B1,B2,B3; indexor',plotangle',EEI'];     
                                    
            end
            
        hell = hell +1;
    
        
        
clearvars timeslice

        
        
        
        
        end
            
    end
       
hell = 1;
length(Array);
no_samples_slice*ntwin;

[ tmp ind] = sort(Array(:,1));Array = Array(ind,:);



fwrite(fid5,Array(:,2),'float32');
fwrite(fid4,Array(:,3),'float32');
clearvars Array


end

zerodata = zeros(no_samples_slice*remainder,1)



fwrite(fid5,zerodata,'float32');
fwrite(fid4,zerodata,'float32');

close(writerObj);

end