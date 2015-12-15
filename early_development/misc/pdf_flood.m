function [] = pdf_flood(file_in,sc)
% Replaces zeros in OpendTect or Rokdoc 2D PDFs with the number sc (small
% constant). Writes the resulting PDF in same format as input file. Output
% filename has "_flooded" appended to whatever the input filename was.

    % Read the headers and PDF data from the file 'file_in'
    header = 13;
    hdrtext = textscan(fopen(file_in),'%s%s','delimiter',':');
    hdrvalue = str2double(hdrtext{2});
    while isnan(hdrvalue(end,1))
        hdrvalue = hdrvalue(1:end-1,1);
    end
    fclose all;
    [pdf] = cell2mat(textscan(fopen(file_in),'','headerlines',header));
    fclose all;  
    
    % Add the small constant 'sc' to all zero values of the PDF
    pdf(pdf==0) = sc;

    % Write the old headers and new PDF data to the file 'file_in'_flooded
    fid = fopen(sprintf('%s_flooded',file_in),'wt');
    fprintf(fid,'Probability Density Function 2D "%s" flood by %f output by Jonathan Edgar Matlab code\n\n',file_in,sc);
    fprintf(fid,'X Log Type        : %s\n',cell2mat(hdrtext{2}(length(hdrvalue)-7)));
    fprintf(fid,'X Mid Bin Minimum : %f\n',hdrvalue(end-6,1));
    fprintf(fid,'X Bin Width       : %f\n',hdrvalue(end-5,1));
    fprintf(fid,'X No of Bins      : %d\n\n',hdrvalue(end-4,1));
    fprintf(fid,'Y Log Type        : %s\n',cell2mat(hdrtext{2}(length(hdrvalue)-3)));
    fprintf(fid,'Y Mid Bin Minimum : %f\n',hdrvalue(end-2,1));
    fprintf(fid,'Y Bin Width       : %f\n',hdrvalue(end-1,1));
    fprintf(fid,'Y No of Bins      : %d\n\n',hdrvalue(end,1));
    [nx ny]=size(pdf);
    for y = 0:ny-1
        fprintf(fid,'X Bin %d\t',y);
    end
    fprintf(fid,'\n');
    for x = 1:nx
        for y = 1:ny
            fprintf(fid,'%f\t',pdf(x,y));
        end
        fprintf(fid,'\n');
    end
    fclose(fid);

end

