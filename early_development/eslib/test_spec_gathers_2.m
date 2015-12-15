clear all
close all

addpath('/apps/gsc/matlab-mcode');
addpath('/TZA/dtect/TZA_Mafia/spec_decomp');
cd '/data/TZA/dtect/TZA_Mafia/spec_decomp'

file_in = 'pre_stack_mcmp_4x4';

file=fopen(file_in);
    [header]=textscan(file,'%d %d %d %*[^\n]',1,'bufsize',8191);
    [geometry]=textscan(file,'%d %d %d %*[^\n]','HeaderLines',1,'bufsize',8191);
fclose(file);

n_offsets = length(unique(geometry{3}));
n_ilines = length(unique(geometry{1}));
n_xlines = length(unique(geometry{2}));
n_samples = header{3};
s_rate = header{2};

dlmwrite('pre_gathers_med_2',header,'delimiter','\t','-append');

%tStart=tic;
for gather_n=1:n_ilines*n_xlines;
    fprintf('Processing gather: %d of %d\n',gather_n,n_ilines*n_xlines);
    skip=(gather_n-1)*n_offsets;
    file=fopen(file_in);
        data_cellarray = textscan(file,'%d',n_offsets*(n_samples+3),'HeaderLines',1+skip,'bufsize',8191);
    fclose(file);

    data = double(reshape(cell2mat(data_cellarray),n_samples+3,n_offsets));
    header_label = data(1:3,:);
    data = data(4:length(data),:);

    [frequency, maximum, gradient, index] = spec_gathers(data, header_label(3,:), 0:s_rate:(n_samples-1)*s_rate);
    
    max_freq = medfilt1(frequency/10,21);
    
    processed_gathers(1,gather_n) = header_label(1,1);
    processed_gathers(2,gather_n) = header_label(2,1);
    processed_gathers(3:n_samples+2,gather_n) = max_freq;    
  
    % dlmwrite('processed_gathers',[processed_gathers','delimiter','\t','-append'); 

    % dlmwrite('pre_gathers_med',[header_label(1,1) header_label(2,1) max_freq'],'delimiter','\t','-append');
end

dlmwrite('pre_gathers_med_2',processed_gathers','delimiter','\t','-append');

% figure(1)
%     plot(frequency,0:s_rate:(n_samples-1)*s_rate)
%     set(gca,'YDir','reverse')

%figure(1)
%imagesc(processed_gathers(2,1):2:processed_gathers(2,n_ilines),0:s_rate:(n_samples-1)*s_rate,processed_gathers(3:n_samples+2,1:inline))
    
 
% plot(frequency/10,0:s_rate:(n_samples-1)*s_rate); set(gca,'YDir','reverse')

%figure(1)
%pcolor(data); %hold on; stem(index,0:4:12); 
%set(gca,'XTick',10:10:40);
%set(gca,'YTick',0:4:12)
