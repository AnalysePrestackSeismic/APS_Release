%segy = segy_read_files('/data/NOR/dtect/gullris_3d_site_survey/Misc');
segy = segy_read_files('C:\Users\James\Work\data');

%get "crossline"
% xline=1;
% xltraces=zeros(1116,573);
% for i=1:533
%     traces=segy_read_traces(segy{1},573*(i-1)+xline,573*(i-1)+xline,0,0);
%     xltraces(:,i)=traces.data;
% end
% figure('Color',[1 1 1]);colorbar;
% contourf(flipud(xltraces));

traces = segy_read_traces(segy{1},1,305409,0,0);
%figure('Color',[1 1 1]);colorbar;
%contourf(flipud(traces.data));
