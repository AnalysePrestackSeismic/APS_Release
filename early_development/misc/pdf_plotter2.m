% matlab code for plotting, smoothing and saving 3D surface and 2D contour pdfs (from rocdoc or opendtect)
% author: Jonathan Edgar

%reset the work area
clear all
close all
clc

%turn off text interpretation for plots (to ensure underscores are not interpreted as latex subscript)
set(0,'DefaulttextInterpreter','none');

%set number of pdfs to plot
numplot = input('Enter number of pdfs to plot (5 maximum): ');

%set single or multiple plots
if numplot > 1
    oneplot = input('Enter 1 to plot all pdfs on the same graph or 0 to plot on seperate graphs: ');
else
    oneplot = 0;
end

if oneplot > 0
    norm = input('Enter 1 to normalise all pdfs to unity (for visualisation only) or 0 maintain true z relationships: ');
else
    norm = 0;
end

%set whether to smooth the pdfs or not
filtornot = input('Enter 1 to apply average smoothers or 0 to not smooth: ');

%define the colours to use for the plots (red, blue, green, yellow, grey)
colours=[1,0,0;0,0,1;0,1,0;1,1,0;0.8,0.8,0.8];

%begin plotting loop
for c=1:numplot
    file_in = input('''Input file name'' = ');                      %select pdf file to plot
    header = input('Number of header lines to skip: ');             %skip header lines (13 if a pdf from rokdoc)
    [hdrtext hdrvalue_tmp] = textread(file_in,'%s%s','delimiter',':');  %read the headers in the file
    hdrvalue = str2double(hdrvalue_tmp);
    while isnan(hdrvalue(end,1))
        hdrvalue = hdrvalue(1:end-1,1);
    end
    hdrvalue_tmp = hdrvalue_tmp(1:length(hdrvalue),1);
    [pdf] = textread(file_in,'','headerlines',header);              %read the pdf matrix in the file

    x_min = hdrvalue(end-6,1);    %set minimum x value defined in the headers
    x_int = hdrvalue(end-5,1);    %set x bin size defined in the headers
    y_min = hdrvalue(end-2,1);    %set minimum y value defined in the headers
    y_int = hdrvalue(end-1,1);    %set y bin size defined in the headers
    
    [nx ny]=size(pdf);                              %calculate size of the pdf
    x_axis = (((nx-1)*x_int)+x_min:-x_int:x_min);   %build the x-axis
    y_axis = (y_min:y_int:((ny-1)*y_int)+y_min);    %build the y-axis
    %x_axis = fliplr(x_axis);
    %y_axis = fliplr(y_axis);
    
    %smoothing filter builder
%     if filtornot > 0
%         filttype=input(sprintf('Select an average smoother to apply to %s (1 = none, 2 = 2x2, 3 = 3x3, etc): ',file_in));   %define filter size
%         if filttype < 2
%             filter=1;                                           %equivalent to not filtering
%         else
%             filter=ones(filttype,filttype)/(filttype*filttype); %build convolutional average filter matrix
%         end
%     else
%         filter=1;                                               %equivalent to not filtering
%         filttype=1;
%     end
% 
%     %make plots
%     if oneplot > 0  %single plot
%         
%         figure(1)
%         for sub=1:3
%             subplot(2,2,sub)
%             if norm > 0
%                 surf(y_axis,x_axis,((conv2(pdf,filter,'same'))/(max(max(conv2(pdf,filter,'same'))))),'FaceAlpha',0.5,'EdgeAlpha',0.5)    %normalised 3D surface plot
%                 title(sprintf('Normalised 3D Surface Plot of PDFs (view %d)',sub))              %title the plot
%                 zlabel('Normalised Number')                                                     %label the z-axis
%             else
%                 surf(y_axis,x_axis,conv2(pdf,filter,'same'),'FaceAlpha',0.5,'EdgeAlpha',0.5)    %3D surface plot
%                 title(sprintf('3D Surface Plot of PDFs (view %d)',sub))                         %title the plot
%                 zlabel('Number')                                                                %label the z-axis
%             end
%             colormap(colours(c,:))                                                              %use colour from colours matrix
%             freezeColors                                                                        %lock the colourmap (requires freezeColors.m)
%             axis tight                                                                          %resize axis
%             ylabel('AI')                                                                        %label the y-axis
%             xlabel('Vp/Vs')                                                                     %label the x-axis
%             if sub==1
%                 view([-37.5 36]);                                                               %rotate the view
%             else
%                 view([(((sub-2)*90)-10) 36]);
%             end
%             hold all                                                                            %hold everything previously plotted
%         end
%                 
%         subplot(2,2,4)
%         contour(y_axis,x_axis,conv2(pdf,filter,'same'),25); %2D contour plot
%         colormap(colours(c,:))                              %use colour from colours matrix
%         freezeColors                                        %lock the colourmap (requires freezeColors.m)
%         axis tight                                          %resize axis
%         ylabel('AI')                                        %label the y-axis
%         xlabel('Vp/Vs')                                     %label the x-axis
%         title('2D Contour Plot of PDFs')                    %title the plot
%         hold all                                            %hold everything previously plotted
%                 
%         if c == numplot                                         %when everything has been plotted save the graph
%             set(gcf, 'PaperUnits', 'inches');                   %set dimensions to inches
%             set(gcf, 'PaperPosition', [1,1,14,12]);             %set picture size
%             print('-r300','-dtiff','overlain_pdf_plots.tiff');  %save as .tiff at 300dpi
%         end
%         
%     else    %multiple plots
%         
%         figure(c)                                                                           %increment figure number
%         for sub=1:3
%             subplot(2,2,sub)
%             surf(y_axis,x_axis,conv2(pdf,filter,'same'),'FaceAlpha',0.8,'EdgeAlpha',0.4)    %3D surface plot
%             axis tight                                                                      %resize axis
%             ylabel('AI')                                                                    %label the y-axis
%             xlabel('Vp/Vs')                                                                 %label the x-axis
%             zlabel('Number')                                                                %label the z-axis
%             title(sprintf('3D Surface Plot of PDFs (view %d)',sub))                         %title the plot
%             if sub==1
%                 view([-37.5 36]);                                                           %rotate the view
%             else
%                 view([(((sub-2)*90)-10) 36]);
%             end
%         end
%         
%         subplot(2,2,4)
%         contourf(y_axis,x_axis,conv2(pdf,filter,'same'),25);    %2D filled contour plot
%         alpha(0.5)                                              %make colour semi-transparent
%         axis tight
%         ylabel('AI')                                            %label the y-axis
%         xlabel('Vp/Vs')                                         %label the x-axis
%         title(sprintf('2D Contour Plot of %s',file_in))         %title the plot with the pdf filename
%         colormap(jet)                                           %use default colours
%         
%         set(gcf, 'PaperUnits', 'inches');                               %set dimensions to inches
%         set(gcf, 'PaperPosition', [1,1,14,12]);                         %set picture size
%         print('-r300','-dtiff',sprintf('plot_%d_%s.tiff',c,file_in));   %save as .tiff at 300dpi
%     end
%     
%     if filttype > 1                                                                         %give option to save a filtered pdf
%         savepdf = input('Enter 1 to save the filtered PDF as a .txt or 0 to not save: ');
%     else
%         savepdf = 0;
%     end
%     
%     if savepdf > 0
%         %code to output .txt file in rokdoc and opendtect format
%         fid = fopen(sprintf('%dx%d_filtered_%s',filttype,filttype,file_in),'wt');
%         fprintf(fid,'Probability Density Function 2D "%dx%d_filtered_%s" output by Jonathan Edgar Matlab code\n\n',filttype,filttype,file_in);
%         fprintf(fid,'X Log Type        : %s\n',hdrvalue{2});
%         fprintf(fid,'X Mid Bin Minimum : %s\n',hdrvalue{3});
%         fprintf(fid,'X Bin Width       : %s\n',hdrvalue{4});
%         fprintf(fid,'X No of Bins      : %s\n\n',hdrvalue{5});
%         fprintf(fid,'Y Log Type        : %s\n',hdrvalue{6});
%         fprintf(fid,'Y Mid Bin Minimum : %s\n',hdrvalue{7});
%         fprintf(fid,'Y Bin Width       : %s\n',hdrvalue{8});
%         fprintf(fid,'Y No of Bins      : %s\n\n',hdrvalue{9});
%         for y = 0:ny-2
%             fprintf(fid,'X Bin %d\t',y);
%         end
%         fprintf(fid,'\n');
%         filt_pdf = conv2(pdf,filter,'same');
%         for x = 1:nx
%             for y = 1:ny
%                 fprintf(fid,'%f\t',filt_pdf(x,y));
%             end
%             fprintf(fid,'\n');
%         end
%         fclose(fid);
%     end
    
    savepdfpetrel = input('Enter 1 to save the PDF as a .txt (Petrel format) or 0 to not save: ');
    
    if savepdfpetrel > 0
        
%         fid1 = fopen(sprintf('column_format_1_%s',file_in),'wt');
%         fprintf(fid1,'Probability Density Function 2D output by Jonathan Edgar Matlab code\n\n');
%         fprintf(fid1,'%s\t\t\t',hdrvalue_tmp{end-7});
%         fprintf(fid1,'%s\t\t',hdrvalue_tmp{end-3});
%         fprintf(fid1,'Probability\n');
%         for x = 1:nx
%             for y = 1:ny
%                 fprintf(fid1,'%f\t%f\t%f\n',x_min+(x-1)*x_int,y_min+(y-1)*y_int,pdf(x,y));
%             end
%         end
%         fclose(fid1);
%         
%         fid2 = fopen(sprintf('column_format_2_%s',file_in),'wt');
%         fprintf(fid2,'Probability Density Function 2D output by Jonathan Edgar Matlab code\n\n');
%         fprintf(fid2,'%s\t\t\t',hdrvalue_tmp{end-7});
%         fprintf(fid2,'%s\t\t',hdrvalue_tmp{end-3});
%         fprintf(fid2,'Probability\n');
%         for y = 1:ny
%             for x = 1:nx
%                 fprintf(fid2,'%f\t%f\t%f\n',x_min+(x-1)*x_int,y_min+(y-1)*y_int,pdf(x,y));
%             end
%         end
%         fclose(fid2);
        
        fid3 = fopen(sprintf('column_format_3_%s',file_in),'wt');
        fprintf(fid3,'Probability Density Function 2D output by Jonathan Edgar Matlab code\n\n');
        fprintf(fid3,'%s\t\t\t',hdrvalue_tmp{end-3});
        fprintf(fid3,'%s\t\t',hdrvalue_tmp{end-7});
        fprintf(fid3,'Probability\n');
        for x = 1:nx
            for y = 1:ny
                fprintf(fid3,'%f\t%f\t%f\n',1000*y_axis(y),x_axis(x),pdf(x,y));
            end
        end
        fclose(fid3);
        
%         fid4 = fopen(sprintf('column_format_4_%s',file_in),'wt');
%         fprintf(fid4,'Probability Density Function 2D output by Jonathan Edgar Matlab code\n\n');
%         fprintf(fid4,'%s\t\t\t',hdrvalue_tmp{end-3});
%         fprintf(fid4,'%s\t\t',hdrvalue_tmp{end-7});
%         fprintf(fid4,'Probability\n');
%         for y = 1:ny
%             for x = 1:nx
%                 fprintf(fid4,'%f\t%f\t%f\n',y_min+(y-1)*y_int,x_min+(x-1)*x_int,pdf(x,y));
%             end
%         end
%         fclose(fid4);
        
    end
    
end