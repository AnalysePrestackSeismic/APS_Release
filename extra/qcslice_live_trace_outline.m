function [live_trace_polygon] = qcslice_live_trace_outline(filename)
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.
%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
% qcslice_live_trace_outline(filename) 
% this reads a xy polygon and zeros out all the data that is not inside
% that polygon and then writes it back out to a new ieee segy file
%%-------------------------------------------------------------------------
% constants

%load the file
load(filename);

%     % Inline / crossline reading step
%     trace_ilxl_bytes(:,1) = int64(trace_header(bytes_to_samples == xloc_byte,:))';
%     trace_ilxl_bytes(:,2) = int64(trace_header(bytes_to_samples == yloc_byte,:))';
%     %
%     % make a couple of matlab arrays to hold a couple of time slices
%     qcslice(1:4,(((ii)*blocktr)+1):end) = [single(trace_ilxl_bytes(:,1)');single(trace_ilxl_bytes(:,2)');traces(floor((60+strseismic.n_samples)/3),:);traces(floor((60+strseismic.n_samples)/2),:)];
%     out_file_name_slice_tmp = strcat(out_file_name_slice,'_qc_out.mat');
%     save(out_file_name_slice_tmp,'qcslice','endcj','-v7.3'); % Saves Seismic structure to mat file
%     
    
    hor_file = [qcslice(1,:)',qcslice(2,:)'];
    
    xorign = 730848.24;
    yorign = 5923958.23;
    xlineinc = 15;
    inlineinc = 12.5;
    firstinline = 8;
    firstxline = 8;
    rotang = 360-41;
    ilnoinc = 1;
    xlnoinc = 1;
    coscal = 0.01;
    
    
    easting = qcslice(1,:)'*coscal;
    northing = qcslice(2,:)'*coscal;
%     easting = [1036214.91; 872265.63; 894797.52];
%     northing = [5989627.60; 5801025.79; 6112560.04];
    
    xdist = (((easting-xorign)*cosd(rotang))+((northing-yorign)*sind(rotang)));
    ydist = (((northing-yorign)*cosd(rotang))-((easting-xorign)*sind(rotang)));
    
    inline = firstinline+(floor((ydist/inlineinc)+0.5)*ilnoinc);
    xline = firstxline+(floor((xdist/xlineinc)+0.5)*xlnoinc);

    hor_file = [inline,xline];
    
    int_n_traces =  size(qcslice,2);

    min_il = min(hor_file(:,1));
    max_il = max(hor_file(:,1));
    min_xl = min(hor_file(:,2));
    max_xl = max(hor_file(:,2));

    il_inc = mode(diff(unique(hor_file(:,1)))); % Primary Key (inline) increment (mode )
    xl_inc = mode(diff(unique(hor_file(:,2)))); % Secondary Key (inline) Increment (mode)

    pkeyn = 1+(max_il-min_il)/il_inc;           % Calculate Number of inlines (primary key)
    skeyn = 1+(max_xl-min_xl)/xl_inc;

    n_iline = (hor_file(:,1)-min_il)/il_inc+1;
    n_xline = (hor_file(:,2)-min_xl)/xl_inc+1;
    lin_ind = ((n_iline-1).*skeyn)+n_xline;

    slice = zeros(skeyn,pkeyn);
    slice(lin_ind) = ones(int_n_traces,1);    
    slice = imfill(slice);    
    slice_x = zeros(skeyn,pkeyn);
    slice_y = zeros(skeyn,pkeyn);
    slice_x(lin_ind) = xline;
    slice_y(lin_ind) = inline;
    
    [B,L,N,A] = bwboundaries(slice,'noholes');
       
    boundary = B{1};    
    
    hor_file = [boundary(:,2),boundary(:,1)];

    min_il = min(hor_file(:,1));
    max_il = max(hor_file(:,1));
    min_xl = min(hor_file(:,2));
    max_xl = max(hor_file(:,2));

    il_inc = mode(diff(unique(hor_file(:,1)))); % Primary Key (inline) increment (mode )
    xl_inc = mode(diff(unique(hor_file(:,2)))); % Secondary Key (inline) Increment (mode)

    pkeyn = 1+(max_il-min_il)/il_inc;           % Calculate Number of inlines (primary key)
    skeyn = 1+(max_xl-min_xl)/xl_inc;

    n_iline = (hor_file(:,1)-min_il)/il_inc+1;
    n_xline = (hor_file(:,2)-min_xl)/xl_inc+1;
    lin_ind_bound = ((n_iline-1).*skeyn)+n_xline;    
    
    live_trace_polygon = [slice_x(lin_ind_bound),slice_y(lin_ind_bound)];
    
    inlines_present = unique(slice_x(lin_ind_bound));
    out_inlines = zeros(size(inlines_present,1),2);
    
    for mm = 1 : size(inlines_present,1)
        no_of_traces_tmp = sort(live_trace_polygon(live_trace_polygon(:,1) == inlines_present(mm),2)');
        diffout = diff(unique(no_of_traces_tmp));
        diffout(diffout > 1) = diffout( diffout > 1)+1;
        if sum(diffout ~= 1) == 0
            out_inlines(mm,:) = [inlines_present(mm), (sum(diffout) +1)];
        else
            out_inlines(mm,:) = [inlines_present(mm), (sum(diffout)+sum( find(diffout ~= 1) > 1))];
        end
    end

    %figure; plot(live_trace_polygon(:,1),live_trace_polygon(:,2),'-r')
    
    filenameouttmp = strsplit(filename,'.');
    filenameout = strcat(char(filenameouttmp(1)),'.ntrc');
    dlmwrite(filenameout,out_inlines);
    
    
end




