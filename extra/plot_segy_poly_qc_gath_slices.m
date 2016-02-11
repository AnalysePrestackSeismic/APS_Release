function plot_segy_poly_qc_gath_slices( filepathin )
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

%PLOT_SEGY_POLY_QC_SLICES


% xv = [
%  962647
%  923218
%  851819
%  905189
%  962647
%  ];
%
% yv = [
%  5929028
%  5883670
%  5945736
%  6007131
%  5929028
%  ];
%
% xv = xv.*100;
% yv = yv.*100;
% xv = int64(xv);
% yv = int64(yv);

pltad = 100;

%
% [stat, gatherfiles] = system('ls /data/URY/segy/2013_pgs_uruguay_processing/full_area_final_deliverables_phase1_and_2/final_full_area_gathers_links/gather_final_presdm_angle_gaths_block*.sgy')
%
% while lpi <= size(gatherfiles{1},1)
%     char(gatherfiles{:}(lpi))
%     lpi = lpi + 1;
% end
% scan the directory for all the files that match the pattern
filename_string = 'slice_qc_out.mat';
[filter_files,nfiles] = directory_scan(filepathin,filename_string);

% make a guess at the fold of the gathers, load first gather qc file and then
% findnumber of dupilcated x and y locations
load([filter_files.path{1} filter_files.names{1}]);
%gathfold = unique(qcslice(2,1:500));
pltd = 1000;
%pltd = pltad *  sum(ismember(qcslice(2,1:500),gathfold(1)));
% [qcsliceidx,ia,ic] = unique(A) 
% qcsliceidx = unique(qcslice(1,1:endcj))
% test with just a few blocks
%blocks_tofind = [238 181 182 44 45 6 157];

% the endcj variable comes from the mat file
% the 6 entries in qcslice are x, y , sample value 1 without the polygon applied, sample value 2 without the polygon applied, sample
% value 1 with polygon restriction applied, sample 2 value with polygon restriction applied

qcslice_plot = [qcslice(1,1:pltd:endcj);qcslice(2,1:pltd:endcj);qcslice(3,1:pltd:endcj);qcslice(5,1:pltd:endcj)];
clear qcslice endcj;
%
for ii=2:1:nfiles
    %for kk=1:size(blocks_tofind,2)
        %if strfind(filter_files.names{ii},num2str(blocks_tofind(kk)))
            load([filter_files.path{ii} filter_files.names{ii}]);
            %gathfold = unique(qcslice(2,1:500));
            %pltd = pltad *  sum(ismember(qcslice(2,1:500),gathfold(1)));
            qcslice_plot = horzcat(qcslice_plot,[qcslice(1,1:pltd:endcj);qcslice(2,1:pltd:endcj);qcslice(3,1:pltd:endcj);qcslice(5,1:pltd:endcj)]);
            clear qcslice endcj;
        %end
    %end
end
%
red = [repmat(1,50,1);(1:-0.02:0)'];
blue = [(0:0.02:1)';repmat(1,50,1)];
green = [(0:0.02:1)';(0.98:-0.02:0)'];
redblue = [red green blue];

scrsz = get(0,'ScreenSize');
cjf = figure('Position',[50 scrsz(4) scrsz(3)/1.1 (scrsz(4)/1.1)],'PaperType','A1','Visible','on','PaperPositionMode','auto','PaperOrientation','landscape');
colormap(redblue);

colrngcj = (0.5*max(qcslice_plot(3,:)));
%props = {'LineStyle','none','Marker','o','MarkerEdge','b','MarkerSize',6};
%line([qcslice(1,1:10:endcj),qcslice(1,1:10:endcj)],[qcslice(2,1:10:endcj),qcslice(2,1:10:endcj)],props{:});
subplot(1,2,1)
%scatter(qcslice_plot(1,:),qcslice_plot(2,:),20,qcslice_plot(3,:),'filled');
dropfac = floor(size(qcslice_plot(1,:),2)/500000)
scatter(qcslice_plot(1,1:dropfac:end),qcslice_plot(2,1:dropfac:end),20,qcslice_plot(3,1:dropfac:end),'filled');
caxis([(-1*colrngcj) colrngcj]);
colorbar;
% hold on
% plot(xv,yv)
% hold off


%out_file_name_slice_tmp = strcat(out_file_name_slice,'_whole.pdf');

%print(cjf,'-dpdf','-noui','-r200',out_file_name_slice_tmp);


%tr_use = inpolygon(int64(qcslice(1,1:pltd:endcj)),int64(qcslice(2,1:pltd:endcj)),xv,yv);
%     %use logical to set the third dimension of the trace length to zero
%     tracesout = bsxfun(@times,traces,tr_use');

subplot(1,2,2)
%scatter(qcslice_plot(1,:),qcslice_plot(2,:),20,qcslice_plot(4,:),'filled');

scatter(qcslice_plot(1,1:dropfac:end),qcslice_plot(2,1:dropfac:end),20,qcslice_plot(4,1:dropfac:end),'filled');

caxis([(-1*colrngcj) colrngcj]);
colorbar;
% hold on
% plot(xv,yv)
% hold off
%colorbar('peer',axes1);
%out_file_name_slice_tmp = strcat(out_file_name_slice,'_sub_polygon.pdf');
%print(cjf,'-dpdf','-noui','-r200',out_file_name_slice_tmp);

%close cjf;
end

