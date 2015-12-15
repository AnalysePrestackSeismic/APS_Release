clear all
close all

addpath('/apps/gsc/matlab-mcode');
addpath('/TZA/dtect/TZA_Mafia/spec_decomp');
cd '/data/TZA/dtect/TZA_Mafia/spec_decomp'

file_in = 'pre_stack_mcmp_4x4';

file=fopen(file_in);
    [header]=textscan(file,'%d %d %d %*[^\n]',1,'bufsize',20191);
    % [geometry]=textscan(file,'%d %d %d %*[^\n]','HeaderLines',1,'bufsize',20191);
fclose(file);

n_offsets = 15; %length(unique(geometry{3}));
n_ilines = 400; %length(unique(geometry{1}));
n_xlines = 2052; %length(unique(geometry{2}));
n_samples = header{3};
s_rate = header{2};

% dlmwrite('pre_gathers_full_survey2',header,'delimiter','\t','-append');

%tStart=tic;
% for gather_n=1:n_ilines*n_xlines;

n_gathers = n_ilines*n_xlines;
gather_step = 5700;
gather_n = 1;

for gather_grp=1:gather_step:11400 %n_gathers;
    fprintf('Processing gather group: %d of %d\n',gather_grp,n_gathers/gather_step);
    while gather_n < gather_grp;
         fprintf('Processing gather: %d of %d\n',gather_n,n_ilines*n_xlines);
         skip=(gather_n-1)*n_offsets;
         file=fopen(file_in);
             data_cellarray = textscan(file,'%d',n_offsets*(n_samples+3),'HeaderLines',1+skip,'bufsize',20191);
         fclose(file);
 
         data = double(reshape(cell2mat(data_cellarray),n_samples+3,n_offsets));
         header_label = data(1:3,:);
         data = data(4:length(data),:);
 
        % [frequency, maximum, gradient, index] = spec_gathers(data, 5:5:75, 0:s_rate:(n_samples-1)*s_rate);
        
        freq_range = 5:5:75;
        time_range = 0:s_rate:(n_samples-1)*s_rate;        
        
        [samples,frequencies] = size(data);

        [maximum, index] = max(data,[],2);

        frequency = freq_range(index)';

        %for i=1:1:samples-1
        %    gradient(i) = (frequency(i+1)-frequency(i))/time_range(3);
        %end

        %gradient = gradient';

         % processed_gathers(1,gather_n) = header_label(1,1);
         % processed_gathers(2,gather_n) = header_label(2,1);
         processed_gathers(1:n_samples,gather_n) = frequency/10;    
         
         % dlmwrite('processed_gathers',[processed_gathers','delimiter','\t','-append'); 

        gather_n = gather_n + 1;
   end
        %max_freq(:,gather_n) = medfilt1(frequency/10,21);

        % dlmwrite('pre_gathers_full_survey_nogeom',[processed_gathers'],'delimiter','\t','-append');
        
        clear frequency
        clear maximum
        clear gradient
        clear index
        clear processed_gathers        
end
    
% function [frequency, maximum, gradient, index] = spec_gathers(gather, freq_range, time_range)

% function to pick peaks of spectral decomposition gathers
% Input: 
%       - a pre-stack spectral decomposition gather
%       - column is frequency
% Output:
%       -
%       - 
%       - change in frequency with time

% NOTE - might need to add time step

% [samples,frequencies] = size(gather);
% 
% [maximum, index] = max(gather,[],2);
% 
% frequency = freq_range(index)';
% 
% for i=1:1:samples-1
%     gradient(i) = (frequency(i+1)-frequency(i))/time_range(3);
% end
% 
% gradient = gradient';



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
