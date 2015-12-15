
file_path = '/media/BGADV034/';

seismic = segy_read_files_old(file_path);
 
nfiles = length(seismic);

for ii = 1:1:nfiles
    traces{ii} = segy_read_traces_old(seismic{ii},1,seismic{ii}.n_traces,0,0);
end

% Plot segy
figure(1)
for ii = 1:1:nfiles 
    subplot(ceil(sqrt(nfiles)),ceil(sqrt(nfiles)),ii); imagesc(traces{ii}.data);
    axis([traces{ii}.pos(2,1) traces{ii}.pos(2,seismic{ii}.n_traces) 0 seismic{ii}.n_samples*seismic{ii}.s_rate/1000])
    xlabel('Crossline');
    ylabel('Time');
    title(['Angle stack ', num2str(seismic{ii}.angle),' degrees']);
end 

% Plot
sample_no = 100;
trace_number = 10;

for ii = 1:1:nfiles 
    amp(ii) = traces{ii}.data(sample_no,trace_number);
    angle(ii) = seismic{ii}.angle;
end 

figure(2)
plot(angle,amp)