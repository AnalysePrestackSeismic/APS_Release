function[] = crossplotter2(Xaxisfile,Yaxisfile,Time_Volume, nt,ntwin,bytespersample, divisions, run_number)

fprintf('\n Reading in Data..\n');
%% READ IN DATA AND DETERMINE FILE SIZE/ NO. OF TRACES

fid1 = fopen(Xaxisfile);
fid2 = fopen(Yaxisfile);
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
sampleskip = floor(nt/divisions)-1
samplesperslice = bytespersample*ntrace1



fid4 = fopen(sprintf('%d_EEI_data',run_number),'a');
fid5 = fopen(sprintf('%d_flat_bk_trend_angle',run_number),'a');

% x = -200:0.1:200;
% 

screen_size = get(0, 'ScreenSize');
% 
writerObj=VideoWriter('PCA_crossplotmovie');
writerObj.FrameRate=10;
open(writerObj);



%% Calculate the number of timesamples


nloops = floor(nt/ntwin);



for k = 1:nloops
 if ftell(fid1) < filesize1     
   
    jump = sampleskip*samplesperslice;
      
  for i = 1:ntwin
        
    timeslice(k).intdata((i-1)*ntrace1+1:i*ntrace1,1)  = fread(fid1,ntrace1,sprintf('%d*float32',ntrace1),jump);
    timeslice(k).graddata((i-1)*ntrace1+1:i*ntrace1,1) = fread(fid2,ntrace2,sprintf('%d*float32',ntrace1),jump);
    timeslice(k).timedata((i-1)*ntrace1+1:i*ntrace1,1) = fread(fid3,ntrace3,sprintf('%d*float32',ntrace1),jump);

   end
    
    %% Set wheeler data to NaN where required   
    
%     for j = 1:ntwin*ntrace1
%     if timeslice(k).intdata(j,1) > 1e29
%         timeslice(k).intdata(j,1) =  NaN;
%         timeslice(k).graddata(j,1) = NaN;
%  
%     end
%     if timeslice(k).graddata(j,1) > 1e29
%         timeslice(k).graddata(j,1) = NaN;
%         timeslice(k).intdata(j,1) =  NaN;
%     end
%     end      
       
    timeslice(k).intdata1  = timeslice(k).intdata(isfinite(timeslice(k).intdata));
    timeslice(k).graddata1 = timeslice(k).graddata(isfinite(timeslice(k).graddata));
           
    timeslice(k).zindex = 0.5*ntwin;     
%    timeslice(k).ztime  =  notnumbermean(timeslice(k).timedata((timeslice(k).zindex-1)*ntrace1+1:timeslice(k).zindex*ntrace1,1));
         
    crossplot(k).chidata  = atan(timeslice(k).intdata./timeslice(k).graddata);
    crossplot(k).backgroundangle = notnumbermean(crossplot(k).chidata);
    crossplot(k).gradient = -1/(tan(crossplot(k).backgroundangle));
    crossplot(k).y = crossplot(k).gradient * timeslice(k).intdata;
    
    PCA(k).X = [timeslice(k).intdata1 timeslice(k).graddata1];
    
    [PCA(k).COEFF,PCA(k).SCORE,PCA(k).latent]= princomp1(PCA(k).X);
    
    % Calculate the EEI attribute.
  
    PCA(k).ANGLE = atan(-1*((PCA(k).COEFF(1,1)./PCA(k).COEFF(2,1))));
         
    PCA(k).EEI = timeslice(k).intdata*cos(PCA(k).ANGLE) + timeslice(k).graddata*sin(PCA(k).ANGLE);  
     
    PCA(k).DEG = 180*(PCA(k).ANGLE/pi);
    
    PCA(k).plotangle(1:ntrace1,1) = PCA(k).DEG;
    
    fwrite(fid5,PCA(k).plotangle,'float32');
   
    fwrite(fid4,PCA(k).EEI,'float32');     
    
    time = mean(timeslice(k).timedata);
    
    
    
    plot(timeslice(k).intdata,timeslice(k).graddata,'+');
    axis([-600 600 -600 600])
    title(sprintf(' Timeslice %4.4f ms ', time));
    
        
%     plot(PCA(k).DEG,k,'+');%timeslice(k).ztime,'+')
%     hold on      
    
    
 end
 
end


% %figure
% %plot(timeslice.zindex,PCA.ANGLE)
% x= -1000:0.1:1000;
%  for k = 1:nt
% %  gcf = figure(1);
% %  set(gcf,'Position', [0 0 (screen_size(3)/2) screen_size(4)] );
% %  axis([-500 500 -500 500])
%  k;
%  figure(1);
%  subplot(2,1,1)
%  axis equal square 
%  axis([-1000 1000 -1000 1000]);
%  plot(timeslice(k).intdata,timeslice(k).graddata ,'+');
%  title(sprintf('Time Slice - %d',timeslice(k).zindex));
%  hold on
%  quiver(0,0,50*(PCA(k).COEFF(1,1)),50*(PCA(k).COEFF(2,1)),'r');
%  hold on
%  quiver(0,0,50*(PCA(k).COEFF(1,2)),50*(PCA(k).COEFF(2,2)),'g');
%  hold off
%  subplot(2,1,2)
%  hist(PCA(k).EEI,x);
%  frame=getframe(figure(1));
%  writeVideo(writerObj,frame);
%  end
% figure
% plot(timeslice.zindex,PCA.ANGLE)
 
close(writerObj);


end


