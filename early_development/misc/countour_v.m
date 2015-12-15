input_dir = {'/data/KEN/segy/2012_l10ab_geotrace/final_velocities_time/migration_rms_vel'};

files_in = directory_scan(input_dir);

index_files = input('Enter numbers of segy files to scan (in bracket [], e.g. [1 3 5]): ');
% this will used for all files
il_byte = 189;
xl_byte = 193;
extra_bytes_to_scan = [181 185];
% Do you want to scan all files that you have selected? [1 - Yes, 0 - No]
scan_files = 1;
% processing with aperture?
%aperture = 1;
nfiles = length(index_files);

% scan the input files
input_segy = segy_read_files(scan_files,index_files,files_in,il_byte,xl_byte,extra_bytes_to_scan);

for i_files = 1:1:nfiles
    
    traces{i_files} = segy_read_traces(input_segy{i_files}.filepaths,input_segy{i_files}.n_samples,input_segy{i_files}.trace_pointers);
    
end

%traces_decimate = traces;

% decimate data, should vectorise
factor = 4;
for i_files = 1:1:nfiles
   n_traces =  size(traces{i_files}.data,2);   
   for i_trace=1:1:n_traces
       traces_decimate{i_files}.data(:,i_trace) = decimate(traces{i_files}.data(:,i_trace),factor);
   end  
end

min_v = min(min(traces_decimate{1}.data));
max_v = max(max(traces_decimate{1}.data));

v = linspace(min_v,max_v,250);

% C is a two-row matrix specifying all the contour lines. Each contour line defined in matrix C begins with a column that contains the value of the contour (specified by v and used by clabel), and the number of (x,y) vertices in the contour line. The remaining columns contain the data for the (x,y) pairs. 
% C = [value1 xdata(1) xdata(2) ... xdata(dim1) value2 xdata(1) xdata(2) ... xdata(dim2)...
%    dim1 ydata(1) ydata(2) ... ydata(dim1) dim2 ydata(1) ydata(2) ... ydata(dim2)];

i = 1;
% find isovelocities
for i_cont=1:1:length(v);
    for i_files = 1:1:nfiles
        contour_line{i_files}.cont = contourc(traces_decimate{i_files}.data,[v(i_cont) v(i_cont)]);
        
        % contour(traces_decimate{i_files}.data,[v(i_cont) v(i_cont)]);
        
        % get data from isovel
        n_vals = length(contour_line{i_files}.cont)-1;
        
        if (n_vals > 10)
        
            %traces_decimate{1}.data(round(contour_line{1}.cont(2,2:end)),round(contour_line2{1}.cont(1,2:end)))

            %input_segy{i_files}.trace_pointers(contour_line{i_files}.cont(1,2:end),4)
            %input_segy{i_files}.trace_pointers(contour_line{i_files}.cont(1,2:end),5)
            
            index = logical(abs(ismember(contour_line{i_files}.cont(1,1:end),v(i_cont))-1));
            
            x{i,i_cont} = input_segy{i_files}.trace_pointers(round(contour_line{i_files}.cont(1,index)),4);
            y{i,i_cont} = input_segy{i_files}.trace_pointers(round(contour_line{i_files}.cont(1,index)),5);
            z{i_cont,i} = contour_line{i_files}.cont(2,index);
            % need to contour z too
            idx = sub2ind(size(traces_decimate{i_files}.data),round(contour_line{i_files}.cont(2,index)),round(contour_line{i_files}.cont(1,index)));
            vals{i_cont,i} = traces_decimate{i_files}.data(idx);
        
            i = i + 1;
            clearvars index idx
            
        end
                
    end  
       n = 4;
       x_d =  decimate(cell2mat(x),factor*n);
       y_d =  decimate(cell2mat(y),factor*n);
       z_d =  decimate(cell2mat(z),factor*n);
       vals_d = decimate(cell2mat(vals),factor*n);
       
       
       %F = TriScatteredInterp(x_d,y_d,vals_d);
       %G = TriScatteredInterp(x_d,y_d,z_d');
       
       %v_i{i_cont} = F(xi,yi);
       %z_i{i_cont} = G(xi,yi);
       
        xmax = max(x_d);
        xmin = min(x_d);
        ymax = max(y_d);
        ymin = min(y_d);

        % find isovelocities

        % contour_line{1} = contourc(traces{1}.data,v);

        % make into volume of data and compute isosurface

        n_points = 7500;

        [xi,yi] = meshgrid(linspace(xmin,xmax,n_points),linspace(ymin,ymax,n_points));

        v_i = griddata(x_d,y_d,vals_d,xi,yi);
        z_i = griddata(x_d,y_d,z_d,xi,yi);
        
        meta_out = strcat('/data/KEN/segy/for_2d_interp','/','dtect_hor_',num2str(i_cont));
        
        output = [xi yi z_i v_i];
        
        save(meta_out,'xi','yi','z_i','v_i');
        
        %dlmwrite(meta_out,output,'delimiter', '\t', 'precision', '%.9f')
          
        
        % v_i{i_cont} = griddata(x_d,y_d,vals_d,xi,yi,'v4');
        
       clearvars x y z vals x_d y_d z_d vals_d contour_line v_i z_i output
       % clearvars values contour_line index
       
       i = 1;
        
end

% start = 1;
% for i_files = 1:1:nfiles
%     if i_files == 1
%         len = length(x{i_files,100});
%     else
%         len = length(x_vals) + length(x{i_files,100});
%     end
%     x_vals(start:len,1) = x{i_files,100};
%     y_vals(start:len,1) = y{i_files,100};
%     v_vals(start:len,1) = vals{i_files,100};
%     z_vals(start:len,1) = z{100,i_files};
%     start = length(x_vals)+1;        
% end
% 
% dtect_val = [x_vals;y_vals;z_vals;v_vals];
% 
% figure(1)
% hold
% for i=1:1:20
%     subplot(4,5,i); imagesc(traces{1,i}.data), hold;
%     a = contour(traces_decimate{i}.data,[v(100) v(100)]);
%     scatter(a(1,2:end),a(2,2:end))
%     hold off
% end
% 
% start = 1;
% for i_files = 1:1:nfiles
%     scatter(x{i_files,100},y{i_files,100})
%     hold all
% end
% hold off
% 
% 
% 
% mesh(xi,yi,[1.6e03 1.7e03]), hold
% mesh(xi,yi,[1.6e03 1.7e03]), hold
% for i = 1:1:20
%     plot3(x{i,10},y{i,10},z{i,10},'o')
% end
% hold off


% for i_files = 1:1:nfiles
%     %hold all
%     %scatter(input_segy{i_files}.trace_pointers(:,4),input_segy{i_files}.trace_pointers(:,5));
%     x(start_i:end_i) = input_segy{i_files}.trace_pointers(:,3);
%     y(start_i:end_i) = input_segy{i_files}.trace_pointers(:,4);
%     z(start_i:end_i) = 
%     v(start_i:end_i) = traces{i_files}.data;
%     
% end
%hold off

