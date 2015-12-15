function[] = crossplotter1(Xaxisfile,Yaxisfile,Time_Volume, nt,bytespersample, divisions, run_number)

fprintf('\n Reading in Data..\n');
%% READ IN DATA AND DETERMINE FILE SIZE/ NO. OF TRACES

fid1 = fopen(Xaxisfile);
fid2 = fopen(Yaxisfile);
fid3 = fopen(Time_Volume);

fseek(fid1, 0, 'eof'); 
filesize1 = ftell(fid1) 
frewind(fid1); 

fseek(fid2, 0, 'eof');
filesize2 = ftell(fid2);
frewind(fid2);

fseek(fid3, 0, 'eof');
filesize3 = ftell(fid3);
frewind(fid3);

ntrace1 = filesize1/(nt*bytespersample)
ntrace2 = filesize2/(nt*bytespersample)
ntrace3 = filesize3/(nt*bytespersample)
sampleskip = floor(nt/divisions)-1
samplesperslice = bytespersample*ntrace1



fid4 = fopen(sprintf('%d_EEI_data',run_number),'a');

% x = -200:0.1:200;
% 

screen_size = get(0, 'ScreenSize');
% 
writerObj=VideoWriter('PCA_crossplotmovie');
writerObj.FrameRate=10;
open(writerObj);

%% INITIALISING VARIABLES
% Used to be loop to divisions
for k = 1:nt
 if ftell(fid1) < filesize1     
   
    count = k;
    jump = sampleskip*samplesperslice;
   
    timeslice(k).intdata  = fread(fid1,ntrace1,sprintf('%d*float32',ntrace1),jump);
    
    timeslice(k).graddata = fread(fid2,ntrace2,sprintf('%d*float32',ntrace1),jump);
    
    timeslice(k).timedata = fread(fid3,ntrace3,sprintf('%d*float32',ntrace1),jump);
    
    timeslice(k).zindex   = notnumbermean(timeslice(k).timedata);
       
    crossplot(k).chidata  = atan(timeslice(k).intdata./timeslice(k).graddata);
    crossplot(k).backgroundangle = notnumbermean(crossplot(k).chidata);
    crossplot(k).gradient = -1/(tan(crossplot(k).backgroundangle));
    crossplot(k).y = crossplot(k).gradient * timeslice(k).intdata;
      
        
    PCA(k).X = [timeslice(k).intdata timeslice(k).graddata];
    
    [PCA(k).COEFF,PCA(k).SCORE,PCA(k).latent]= princomp1(PCA(k).X);
    
    % Calculate the EEI attribute.
    %  
  
    PCA(k).ANGLE = atan(-1*((PCA(k).COEFF(1,1)./PCA(k).COEFF(2,1))));
    %     
    PCA(k).EEI = timeslice(k).intdata*cos(PCA(k).ANGLE) + timeslice(k).graddata*sin(PCA(k).ANGLE);
    %        
    fwrite(fid4,PCA(k).EEI,'float32'); 

    PCA(k).DEG = 180*(PCA(k).ANGLE/pi);
   
%     plot(PCA(k).DEG,timeslice(k).zindex,'+')
%     hold on 
% %     


 end
 
end
% %figure
% %plot(timeslice.zindex,PCA.ANGLE)
% 
% x= -1000:0.1:1000;
%  for k = 1:nt
%  
% %  gcf = figure(1);
% %  set(gcf,'Position', [0 0 (screen_size(3)/2) screen_size(4)] );
% %  axis([-500 500 -500 500])
%  k;
% 
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
%  
%  subplot(2,1,2)
%  hist(PCA(k).EEI,x);
%  frame=getframe(figure(1));
%  writeVideo(writerObj,frame);
%   
%  end
% 

% figure
% plot(timeslice.zindex,PCA.ANGLE)
 
close(writerObj);


end


