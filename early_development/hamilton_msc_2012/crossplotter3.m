function[] = crossplotter3(Xaxisfile,Yaxisfile,Time_Volume, nt,ntwin,bytespersample, divisions, run_number, no_xlines,no_bins)

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



fid4 = fopen(sprintf('%d_EEI_data',run_number),'a');
fid5 = fopen(sprintf('%d_flat_bk_trend_angle',run_number),'a');



%% Video Writing Details

screen_size = get(0, 'ScreenSize');
writerObj=VideoWriter('PCA_crossplotmovie');
writerObj.FrameRate=10;
open(writerObj);


%% Do the Crossplotting on binned data rather than whole slices.



nloops = floor(nt/ntwin);


binsize = floor(no_xlines)/(no_bins)

remainder = no_xlines-(no_bins*binsize)


no_samples= (filesize1/4);
no_samples_slice = (filesize1/4)/nt
samples_per_inline = no_xlines

    for k = 1:nt
    
    buffer(:,1) = fread(fid1,no_samples_slice,'float32');
    buffer2(:,1)= fread(fid2,no_samples_slice,'float32');
    count = 1;

    for i = 1:(no_samples_slice/binsize)      


        timeslice(k).graddata(1:binsize,i) =  buffer(((count-1)*binsize)+1:(count*binsize),1);
            timeslice(k).intdata(1:binsize,i)  =  buffer2(((count-1)*binsize)+1:(count*binsize),1);
                count = count + 1;
        
           
%         for j = 1:binsize
%             if  timeslice(k).intdata(j,1) > 1e29
%                     timeslice(k).intdata(j,i) =  NaN;
%                         timeslice(k).graddata(j,i) = NaN;
%  
%             end
%                 if timeslice(k).graddata(j,i) > 1e29
%                      timeslice(k).graddata(j,i) = NaN;
%                         timeslice(k).intdata(j,i) =  NaN;
%                 end
%         end
        


     
        
% timeslice(k).intdata1(1:binsize,i)  = timeslice(k).intdata(isfinite(timeslice(k).intdata(1:binsize,i)));
% timeslice(k).graddata1(1:binsize,i) = timeslice(k).graddata(isfinite(timeslice(k).graddata(1:binsize,i))); 
PCA(k).X = [timeslice(k).intdata(1:binsize,i) timeslice(k).graddata(1:binsize,i)]; 
[PCA(k).COEFF,PCA(k).SCORE,PCA(k).latent]= princomp1(PCA(k).X);  
PCA(k).ANGLE = atan(-1*((PCA(k).COEFF(1,1)./PCA(k).COEFF(2,1))));  
PCA(k).EEI = timeslice(k).intdata(1:binsize,i)*cos(PCA(k).ANGLE) + timeslice(k).graddata(1:binsize,i)*sin(PCA(k).ANGLE);   
PCA(k).DEG(i) = 180*(PCA(k).ANGLE/pi);   
PCA(k).plotangle(1:binsize) = PCA(k).DEG(i);


fwrite(fid5,PCA(k).plotangle,'float32');
fwrite(fid4,PCA(k).EEI,'float32');


% if i ==4
%     
% for l = 1:4
%     
% plot(PCA(k).DEG(1),k,'o');
% plot(PCA(k).DEG(2),k,'r+');
% plot(PCA(k).DEG(3),k,'g+');
% plot(PCA(k).DEG(4),k,'c+');
% 
% %timeslice(k).ztime,'+')
% hold on 
% 
% end
% 
% end

    end
        
  
% 
%    [PCA(k).COEFF,PCA(k).SCORE,PCA(k).latent]= princomp1(PCA(k).X);
% %     
% %  Calculate the EEI attribute.
% %   
%    PCA(k).ANGLE = atan(-1*((PCA(k).COEFF(1,1)./PCA(k).COEFF(2,1))));
% %          
%    PCA(k).EEI = timeslice(k).intdata*cos(PCA(k).ANGLE) + timeslice(k).graddata*sin(PCA(k).ANGLE);  
% %      
%    PCA(k).DEG = 180*(PCA(k).ANGLE/pi);
% %     
%    PCA(k).plotangle(1:ntrace1,1) = PCA(k).DEG;
% %     
%    fwrite(fid5,PCA(k).plotangle,'float32');
% %    
%    fwrite(fid4,PCA(k).EEI,'float32');     
% %     
%    plot(PCA(k).DEG,k,'+');%timeslice(k).ztime,'+')
%    hold on      
%     
      
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


