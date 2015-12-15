function segy_write_traces(results_path)

    start_point = pwd;
    cd(results_path);

    % Figure out the number of files in the current directory
    [~,nfiles] = system('ls -B *_process_files.mat | wc -l'); 
    nfiles = str2double(nfiles);
    % Read all filenames and convert from ascii to double
    [~,fnames] = system('ls -B1 *_process_files.mat');
    numeric = double(fnames);
    % Preallocate memory for some variables
    count = 2;
    fname_index = zeros(1,nfiles+1);
    % Loop to separate out each file name from one long character string
    for ij= 1:length(fnames)
        if numeric(1,ij) == 10
            fname_index(1,count) = ij;
            count = count+1;
        end
    end
    % Loop to read each file as ascii into a cell array
    for ik=1:nfiles
        process_files_names{ik} = fnames(1,(fname_index(1,ik)+1):(fname_index(1,ik+1)-1));
        func_name = regexp(process_files_names{ik}, '_process_files.mat', 'split');
        system_for = sprintf('ls -B %s*.gy',func_name{1});
        [~,n_segy] = system(system_for); 
        n_segy = str2double(n_segy);
        if n_segy == 0 || isnan(n_segy)
            process_files_mat = strcat(func_name{1},'_process_files.mat'); 
            %load(process_files_mat);
            process_files = load(process_files_mat,'-mat');
            system_for = strcat(sprintf('ls -B %s_results_block* | wc -l',func_name{1}));
            [~,n_blocks] = system(system_for);
            n_blocks = str2double(n_blocks);
            if n_blocks == process_files.n_blocks 
                output_segy(results_path,func_name{1});
            else
                fprintf('Please check all results have completed. \n');
            end
        else
            fprintf('Segy already exists. Please remove files if you want to re-create them. \n');
        end
    end
       cd(start_point);   
end

function output_segy(results_path,func_name)
    fprintf('Building headers ...\n');
    % Load process files mat for function results
    process_files_mat = strcat(results_path,func_name,'_process_files.mat');
    %load(process_files_mat,'-mat')
    process_files = load(process_files_mat,'-mat');

    % Create file of correct size
    o_format='ieee-be';
    filename = '.temp2.segy';
    fid=fopen(filename,'w',o_format);

    % Make and write EBCIDIC (3200 bytes)
    ascii_header=make_header(func_name); % take input arguments and populate
    [nh1,mh1]=size(ascii_header);
        if nh1 > mh1
            ascii_header=ascii_header';
        end
    fwrite(fid,ascii2ebcdic(ascii_header),'uchar');

    % Write binary reel header
    two_bytes=zeros(194,1);
    dt = process_files.s_rate;
    n_samples = process_files.n_samples;
    datafmt = 1;	% Data-format code
    two_bytes(3:7)=[dt,dt,n_samples,n_samples,datafmt]';
    fwrite(fid,[10 20 30],'int32');
    fwrite(fid,two_bytes,'int16');
    
    % Retype 4-byte headers to float32
    % These headers can be written as float32, but will give correct values
    % when read as int32
    %float32_trace_headers4 = reshape(float32_trace_headers4,[],1);
    %float32_trace_headers4 = dec2bin(float32_trace_headers4,32);
    %float32_trace_headers4 = reshape(trace_headers4,[],1);
%     float32_trace_headers4 = process_files.processing_grid.ilxl_grid(:,1:2)';
%     float32_trace_headers4 = float32_trace_headers4(:);
%     float32_trace_headers4 = dec2bin(float32_trace_headers4,32);
%     sign = bin2dec(float32_trace_headers4(:,1));
%     exponent = bin2dec(float32_trace_headers4(:,2:9));
%     mantissa = bin2dec(float32_trace_headers4(:,10:end));
%     float32_trace_headers4 = ((-1).^sign).*(2.^(exponent-127)).*(1+(mantissa.*(2^-23)));
%     float32_trace_headers4(and(exponent == 0,mantissa ~= 0)) = ((-1).^sign(and(exponent == 0,mantissa ~= 0))).*(2.^(exponent(and(exponent == 0,mantissa ~= 0))-127)).*(2*(mantissa(and(exponent == 0,mantissa ~= 0)).*(2^-23)));
%     float32_trace_headers4(and(mantissa == 0,exponent == 0)) = 0;
%     %float32_trace_headers4 = reshape(float32_trace_headers4,[],process_files.n_iline*process_files.n_xline);
%     float32_trace_headers4 = reshape(float32_trace_headers4,[],2);

    fid_tmp_trace_headers4 = fopen('.tmp_trace_headers4','w');
    fwrite(fid_tmp_trace_headers4,process_files.processing_grid.ilxl_grid(:,1:2),'int32');
    fclose(fid_tmp_trace_headers4);
    fid_tmp_trace_headers4 = fopen('.tmp_trace_headers4','r');
%     float32_trace_headers4 = fread(fid_tmp_trace_headers4,'float32');
    float32_trace_headers4 = fread(fid_tmp_trace_headers4,'*single');
    fclose(fid_tmp_trace_headers4);
    float32_trace_headers4 = reshape(float32_trace_headers4,[],2);

    % Create 4-byte trace headers
    trace_headers4 = single(zeros(60,process_files.n_iline*process_files.n_xline));
    byte = process_files.ilbyte;
    index4 = ((byte-1)/4)+1;
    %trace_headers4(index4,:) = process_files.processing_grid.ilxl_grid(:,1)';
    trace_headers4(index4,:) = float32_trace_headers4(:,1)';
    byte = process_files.xlbyte;
    index4 = ((byte-1)/4)+1; 
    %trace_headers4(index4,:) = process_files.processing_grid.ilxl_grid(:,2)';
    trace_headers4(index4,:) = float32_trace_headers4(:,2)';
    
    % Create 2-byte trace headers
    %trace_headers2 = zeros(120,process_files.n_iline*process_files.n_xline);
    trace_headers2 = single(zeros(120,1));
    byte = 115; % n samples
    index2 = ((byte-1)/2)+1;
    trace_headers2(index2,:) = single(process_files.n_samples);
    byte = 117; % s rate
    index2 = ((byte-1)/2)+1;
    trace_headers2(index2,:) = single(process_files.s_rate);
    
    % Retype 2-byte headers to float32
    % These headers can be written as float32, but will give correct values
    % when read as int16
    %float32_trace_headers2 = reshape(trace_headers2,[],1);
    float32_trace_headers2 = dec2bin(trace_headers2,16);
    float32_trace_headers2_oddidx = logical(reshape([ones(1,60);zeros(1,60)],[],1));
    float32_trace_headers2_evenidx = logical(reshape([zeros(1,60);ones(1,60)],[],1));
    float32_trace_headers2 = strcat(float32_trace_headers2(float32_trace_headers2_evenidx,:),float32_trace_headers2(float32_trace_headers2_oddidx,:));
    sign = bin2dec(float32_trace_headers2(:,1));
    exponent = bin2dec(float32_trace_headers2(:,2:9));
    mantissa = bin2dec(float32_trace_headers2(:,10:end));
    float32_trace_headers2 = ((-1).^sign).*(2.^(exponent-127)).*(1+(mantissa.*(2^-23)));
    float32_trace_headers2(and(exponent == 0,mantissa ~= 0)) = ((-1).^sign(and(exponent == 0,mantissa ~= 0))).*(2.^(exponent(and(exponent == 0,mantissa ~= 0))-127)).*(2*(mantissa(and(exponent == 0,mantissa ~= 0)).*(2^-23)));
    float32_trace_headers2(and(mantissa == 0,exponent == 0)) = 0;
    %float32_trace_headers2 = reshape(float32_trace_headers2,[],process_files.n_iline*process_files.n_xline);
    float32_trace_headers2 = repmat(float32_trace_headers2(float32_trace_headers2~=0),1,process_files.n_iline*process_files.n_xline);
    
    trace_headers2 = single(zeros(60,process_files.n_iline*process_files.n_xline));
    byte = 115; % n samples
    index2 = floor(((byte-1)/4)+1);
    trace_headers2(index2,:) = float32_trace_headers2(1,:);
    byte = 117; % s rate
    index2 = floor(((byte-1)/4)+1);
    trace_headers2(index2,:) = float32_trace_headers2(2,:);
    
    float32_trace_headers = [trace_headers2(1:45,:);trace_headers4(46:60,:)];
    
    % Create file
    % fwrite(fid,[trace_headers4; zeros(process_files.n_samples,process_files.n_iline*process_files.n_xline)],'float32');
    fprintf('Writing dummy segy ...\n');
    tmp_count = 1;
    tmp_n_blocks = 10;
    tmp_step = floor(process_files.n_iline*process_files.n_xline/(tmp_n_blocks-1));
    tmp_leftover = process_files.n_iline*process_files.n_xline - tmp_step*(tmp_n_blocks-1);
    for i_trace = 1:tmp_step:process_files.n_iline*process_files.n_xline
        fprintf('- Writing dummy segy block %d of %d\n',tmp_count,tmp_n_blocks);
%         fwrite(fid,[trace_headers4(:,i_trace); zeros(process_files.n_samples,1)],'float32');
%         fwrite(fid,trace_headers2(1:90,i_trace),'int16'); 
%         fwrite(fid,trace_headers4(46:60,i_trace),'int32');
%         fwrite(fid,zeros(process_files.n_samples,1),'float32');
        fwrite(fid,single(zeros(60+process_files.n_samples,tmp_step)),'single');
        tmp_count = tmp_count+1;
    end
    fprintf('- Writing dummy segy block %d of %d\n',tmp_count,tmp_n_blocks);
    fwrite(fid,single(zeros(60+process_files.n_samples,tmp_leftover)),'single');
    fclose(fid);
    % write headers
%     fid=fopen(filename,'r+',o_format);
%     skip = 3600;
%     fseek(fid,skip,'cof');
%     for i_trace = 1:1:process_files.n_iline*process_files.n_xline
%         fwrite(fid,trace_headers2(1:90,i_trace),'int16'); 
%         fwrite(fid,trace_headers4(46:60,i_trace),'int32');
%         skip = process_files.n_samples;
%         fseek(fid,skip*4,'cof');
%     end
%     fclose(fid);

    switch process_files.distribute_type;
        case 'trace'
            % trace headers change for each block
            results_mat = strcat(results_path,func_name,'_results_block_1.mat');
            load(results_mat,'-mat');

             % Create required segy files
            n_segy = size(results_out,1);

            for i_segy = 1:1:n_segy 
                fprintf('Writing segy for result %s ...\n',results_out{i_segy,1}); 
                for i_block = 1:1:process_files.n_blocks
                    fprintf('- Writing block %d\n',i_block);
                    results_mat = strcat(results_path,func_name,'_results_block_',num2str(i_block),'.mat');
                    load(results_mat,'-mat');

                    positions_mat = strcat(results_path,func_name,'_positions_block_',num2str(i_block),'.mat');
                    process_positions = load(positions_mat,'-mat');

                    if i_block == 1
                        fid=fopen(filename,'r+',o_format);
                        skip = 3600;
                        fseek(fid,skip,'cof');
                    end
                    
                    if i_block == process_files.n_blocks
                        [~, n_col_last] = size(results_out{i_segy,2});
                        fwrite(fid,[float32_trace_headers(:,1+(i_block-1)*n_col:n_col_last+(i_block-1)*n_col);single(results_out{i_segy,2})],'single');
                    else
                        [~, n_col] = size(results_out{i_segy,2});
                        fwrite(fid,[float32_trace_headers(:,1+(i_block-1)*n_col:i_block*n_col);single(results_out{i_segy,2})],'single');
                    end

%                     for i_trace = 1:1:size(process_positions.ilxl_grid,1)
%                         skip = 240;
%                         fseek(fid,skip,'cof');
%                         fwrite(fid,double(num2ibm(results_out{i_segy,2}(:,i_trace))),'uint32'); 
%                     end
                end

                if i_segy == n_segy
                    fclose(fid); 
                    system_mv = sprintf('mv %s %s_%s.segy',filename,func_name,results_out{i_segy,1}); 
                    system(system_mv); 
                else
                    fclose(fid);
                    system_cp = sprintf('cp %s %s_%s.segy',filename,func_name,results_out{i_segy,1}); 
                    system(system_cp); 
                end
            end
        case 'either' 
            % trace headers are same for all blocks

            results_mat = strcat(results_path,func_name,'_results_block_1.mat');
            load(results_mat,'-mat');

             % Create required segy files
             n_segy = size(results_out,1);

             for i_segy = 1:1:n_segy 
                fprintf('-- Writing segy for result %s --\n',results_out{i_segy,1}); 
                for i_block = 1:1:process_files.n_blocks
                    fprintf('- Writing block %d\n',i_block);
                    results_mat = strcat(results_path,func_name,'_results_block_',num2str(i_block),'.mat');
                    load(results_mat,'-mat');

                    positions_mat = strcat(results_path,func_name,'_positions_block_',num2str(i_block),'.mat');
                    process_positions = load(positions_mat,'-mat');

                    fid = fopen(filename,'r+',o_format);
                    if i_block == 1
                        skip = 3600 + 240;
                        fseek(fid,skip,'cof');
                        prev_nz = 0;
                    else
                        skip = 3600+240+prev_nz*4;
                        fseek(fid,skip,'cof');
                    end
                    
                    nz = size(process_positions.z_grid,1);
                    for i_trace = 1:1:process_files.n_iline*process_files.n_xline
                        fwrite(fid,double(num2ibm(results_out{i_segy,2}(1:nz,i_trace))),'uint32');
                        skip = 240+(process_files.n_samples-nz)*4;
                        fseek(fid,skip,'cof');
                    end
                    prev_nz = prev_nz + nz;
                    fclose(fid);
                end
             
                if i_segy == n_segy
                   system_mv = sprintf('mv %s %s_%s.segy',filename,func_name,results_out{i_segy,1}); 
                   system(system_mv); 
                else
                    system_cp = sprintf('cp %s %s_%s.segy',filename,func_name,results_out{i_segy,1}); 
                    system(system_cp); 
                end 
             end
         case 'slice' 
            % trace headers are same for all blocks

            results_mat = strcat(results_path,func_name,'_results_block_1.mat');
            load(results_mat,'-mat');

             % Create required segy files
             n_segy = size(results_out,1);

             for i_segy = 1:1:n_segy 
                fprintf('-- Writing segy for result %s --\n',results_out{i_segy,1}); 
                for i_block = 1:1:process_files.n_blocks
                    fprintf('- Writing block %d\n',i_block);
                    results_mat = strcat(results_path,func_name,'_results_block_',num2str(i_block),'.mat');
                    load(results_mat,'-mat');

                    positions_mat = strcat(results_path,func_name,'_positions_block_',num2str(i_block),'.mat');
                    process_positions = load(positions_mat,'-mat');

                    fid = fopen(filename,'r+',o_format);
                    if i_block == 1
                        skip = 3600 + 240;
                        fseek(fid,skip,'cof');
                        prev_nz = 0;
                    else
                        skip = 3600+240+prev_nz*4;
                        fseek(fid,skip,'cof');
                    end
                    
                    nz = size(process_positions.z_grid,1);
                    for i_trace = 1:1:process_files.n_iline*process_files.n_xline
                        fwrite(fid,double(num2ibm(results_out{i_segy,2}(1:nz,i_trace))),'uint32');
                        skip = 240+(process_files.n_samples-nz)*4;
                        fseek(fid,skip,'cof');
                    end
                    prev_nz = prev_nz + nz;
                    fclose(fid);
                end
             
                if i_segy == n_segy
                   system_mv = sprintf('mv %s %s_%s.segy',filename,func_name,results_out{i_segy,1}); 
                   system(system_mv); 
                else
                    system_cp = sprintf('cp %s %s_%s.segy',filename,func_name,results_out{i_segy,1}); 
                    system(system_cp); 
                end 
             end
    end  
end

function ascii_header=make_header(func_name)
% Function creates ASCII version of standard EBCIC header of SEG-Y format
    ascii_header=char(...
    'C 1 CLIENT                        COMPANY                       CREW NO', ...         
    'C 2 LINE            AREA                        MAP ID ', ...                          
    'C 3 REEL NO           DAY-START OF REEL     YEAR      OBSERVER', ...                   
    'C 4 INSTRUMENT: MFG            MODEL            SERIAL NO', ...                       
    'C 5 DATA TRACES/RECORD        AUXILIARY TRACES/RECORD         CDP FOLD', ...           
    'C 6 SAMPLE INTERNAL         SAMPLES/TRACE       BITS/IN      BYTES/SAMPLE', ...        
    'C 7 RECORDING FORMAT        FORMAT THIS REEL        MEASUREMENT SYSTEM', ...           
    'C 8 SAMPLE CODE: FLOATING PT     FIXED PT     FIXED PT-GAIN     CORRELATED ', ...      
    'C 9 GAIN  TYPE: FIXED     BINARY     FLOATING POINT     OTHER ', ...                   
    'C10 FILTERS: ALIAS     HZ  NOTCH     HZ  BAND    -     HZ  SLOPE    -    DB/OCT ', ...  
    'C11 SOURCE: TYPE            NUMBER/POINT        POINT INTERVAL', ...                   
    'C12     PATTERN:                           LENGTH        WIDTH', ...                   
    'C13 SWEEP: START     HZ  END     HZ  LENGTH      MS  CHANNEL NO     TYPE', ...         
    'C14 TAPER: START LENGTH       MS  END LENGTH       MS  TYPE', ...                      
    'C15 SPREAD: OFFSET        MAX DISTANCE        GROUP INTERVAL', ...                    
    'C16 GEOPHONES: PER GROUP     SPACING     FREQUENCY     MFG          MODEL', ...        
    'C17     PATTERN:                           LENGTH        WIDTH', ...                   
    'C18 TRACES SORTED BY: RECORD     CDP     OTHER', ...                                  
    'C19 AMPLITUDE RECOVEY: NONE      SPHERICAL DIV       AGC    OTHER', ...                
    'C20 MAP PROJECTION                      ZONE ID       COORDINATE UNITS', ...           
    'C21 PROCESSING:', ...                                                                  
    'C22 PROCESSING:', ...                                                                  
    'C23 ', ...                                                                             
    'C24 ', ...                                                                             
    'C25 ', ...                                                                             
    'C26 ', ...                                                                             
    'C27 ', ...                                                                             
    'C28 ', ...                                                                             
    'C29 ', ...                                                                             
    'C30 ', ...                                                                             
    'C31 ', ...                                                                             
    'C32 ', ...                                                                             
    'C33 ', ...                                                                             
    'C34 ', ...                                                                             
    'C35 ', ...                                                                             
    'C36 ', ...                                                                             
    'C37 ', ...                                                                             
    'C38 ', ...                                                                             
    'C39 ', ...                                                                             
    'C40 END EBCDIC')';
end

function ebcdic=ascii2ebcdic(ascii)
% Function converts ASCII string to EBCDIC
% see http://www.room42.com/store/computer_center/code_tables.shtml
% Date Feb. 20, 2000;  written by E. Rietsch
% INPUT
% ascii         ASCII string
% OUTPUT
% ebcdic	EBCDIC string
%   		  ebcdic=ascii2ebcdic(ascii)

    pointer =  ...
      [0  16  64 240 124 215 125 151  75  75  75  75  75  75  75  75
       1  17  90 241 193 216 129 152  75  75  75  75  75  75  75  75
       2  18 127 242 194 217 130 153  75  75  75  75  75  75  75  75
       3  19 123 243 195 226 131 162  75  75  75  75  75  75  75  75
       4  20  91 244 196 227 132 163  75  75  75  75  75  75  75  75
       5  21 108 245 197 228 133 164  75  75  75  75  75  75  75  75
       6  22  80 246 198 229 134 165  75  75  75  75  75  75  75  75
       7  23 125 247 199 230 135 166  75  75  75  75  75  75  75  75
       8  24  77 248 200 231 136 167  75  75  75  75  75  75  75  75
       9  25  93 249 201 232 137 168  75  75  75  75  75  75  75  75
      10  26  92 122 209 233 145 169  75  75  75  75  75  75  75  75
      11  27  78  94 210 173 146 192  75  75  75  75  75  75  75  75
      12  28 107  76 211 224 147 106  75  75  75  75  75  75  75  75
      13  29  96 126 212 189 148 208  75  75  75  75  75  75  75  75
      14  30  75 110 213  95 149 161  75  75  75  75  75  75  75  75
      15  31  97 111 214 109 150  75  75  75  75  75  75  75  75  75];
    pointer=pointer(:);

    ebcdic=pointer(ascii+1);
end

function b=num2ibm(x)

% num2ibm : convert IEEE 754 doubles to IBM 32 bit floating point format
%    b=num2ibm(x)
% x is a matrix of doubles
% b is a corresponding matrix of uint32
%
% The representations for NaN and inf are arbitrary
%
% See also ibm2num

% 
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%
% (C) Brian Farrelly, 22 October 2001
%  mailto:Brian.Farrelly@nho.hydro.com          Norsk Hydro Research Centre
%  phone +47 55 99 68 74                 (((                  Postboks 7190
%  fax   +47 55 99 69 70                2oooS                 N-5020 Bergen
%  home  +47 55 13 78 49                HYDRO                        Norway
%

b=repmat(uint32(0),size(x));
% err=zeros(size(x));

%format long

x(x> 7.236998675585915e+75)= inf;        % change big numbers to infinity
x(x<-7.236998675585915e+75)=-inf;        % 7.236998675585915e+75 is
                                         %    ibm2num(uint32(hex2dec('7fffffff')) or
                                         %    ibm2num(num2ibm(inf))

[F E]=log2(abs(x));

e=E/4;                         % exponent of base 16
ec=ceil(e);                    % adjust upwards to integer
p=ec+64;                       % offset exponent

f=F.*2.^(-4*(ec-e));           % correct mantissa for fractional part of exponent
f=round(f*2^24);               % convert to integer. Roundoff here can be as large as
                               % 0.5/2^20 when mantissa is close to 1/16 so that
                               % 3 bits of signifance are lost.

p(f==2^24)=p(f==2^24)+1;       % Roundoff can cause f to be 2^24 for numbers just under a
f(f==2^24)=2^20;               % power of 16, so correct for this

%format hex
psi=uint32(p*2^24);            % put exponent in first byte of psi.
phi=uint32(f);                 % put mantissa into last 3 bytes of phi 

% make bit representation

b=bitor(psi,phi);                        % exponent and mantissa
b(x<0)=bitset(b(x<0),32);                % sign bit 
%format long

% special cases

b(x==0)          =uint32(0)                  ;         %  bias is incorrect for zero 
b(isnan(x))      =uint32(hex2dec('7fffffff'));         %  7.237005145973116e+75 in IBM format
b(isinf(x) & x>0)=uint32(hex2dec('7ffffff0'));         %  7.236998675585915e+75    ,,
b(isinf(x) & x<0)=uint32(hex2dec('fffffff0'));         % -7.236998675585915e+75    ,,
                                                       % Note that NaN > inf in IBM format

% check bit representation for normal cases 

%checkx=ibm2num(b);                      % note that use of base 16 in IBM format
%z=x==0;                                 % can lead to a loss of 3 bits of precision
%err(z)=0;                               % compared with an IEEE single.
%q=(checkx(~z)-x(~z))./x(~z);
%err(~z) = abs(q) > 5e-7;                % this is almost reached with numbers
                                        % of the form 16^n + 0.5*16^(n-5) where
                                        % the mantissa is 100001 hex. Roundoff
                                        % error is then 0.5/16^5=0.5/2^20=4.7684e-7
                                                    
                                          
% if any(err)
%    warning('Conversion error in num2ibm for the following:')
%    disp(x(logical(err)))
% end   

end