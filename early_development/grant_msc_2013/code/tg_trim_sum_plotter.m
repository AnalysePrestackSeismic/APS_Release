function [output, output_2] = trim_sum_plotter(filepath)
%% A Function that takes each ASCII output from Opendtect and runs them in matlab
% The function reads in the ASCII files and calculates a histogram for each
% set of data specified by the user and outputs the modal bin range, Kurtosis and Skewness of each set
% The code assumes that the dataset is loaded in an N x 5 matrix, where
% columns 1 - 5 correspond to the inline number, crossline number, depth,
% starting model data and final model data, respectively 
% Author: Troy Grant
% Date: 08/07/2013

%%
% Changing the current MALTAB directory to the one specified by the
% filepath
filepath = '/data/NOR/dtect/nor_hegg_fwi_depth';
addpath(filepath);

% Loading the matrix containing the inline and crossline ranges for each
% file (required for the file_name string)
load il_xl_data_range.mat;

% Counting the number of file ranges (and hence, files) in the matrix
no_files = size(il_xl_data_range,1);

for ii = 2:no_files;
file_name = fullfile(sprintf('tg_avg_il%d-%d_xl%d-%d_RMS',il_xl_data_range(ii,1),il_xl_data_range(ii,2),il_xl_data_range(ii,3),il_xl_data_range(ii,4)));
% file_name = fullfile(sprintf('tg_avg_RMS'));
file = dlmread(file_name,'\t');

% Remove NANs
log2 = file(:,5) ~= 1.0000e+30;
log1 = file(:,4) ~= 1.0000e+30;

% Separating the starting and final model values (columns 4 & 5,
% respectively) 
tg_avg_RMS_col4 = file(log1,1:4);
tg_avg_RMS_col5 = [file(log2,1:3) file(log2,5)];

option = input('Would you like to use a subsample of the data? (yes = 1, no = 0) ');

if option == 1
% Asking the user to input the inline and crossline range for subsample
    starting_il_range(:,1) = input('Please specify the min inline you wish to sample ');
    starting_il_range(:,2) = input('Please specify the max inline you wish to sample ');
    final_il_range(:,1) = starting_il_range(:,1);
    final_il_range(:,2) = starting_il_range(:,2);
    
    starting_xl_range(:,1) = input('Please specify the min crossline you wish to sample ');
    starting_xl_range(:,2) = input('Please specify the max crossline you wish to sample ');
    final_xl_range(:,1) = starting_xl_range(:,1);
    final_xl_range(:,2) = starting_xl_range(:,2);
    
else
    starting_il_range(:,1) = min(tg_avg_RMS_col4(:,1));
    starting_il_range(:,2) = max(tg_avg_RMS_col4(:,1));
    final_il_range(:,1) = min(tg_avg_RMS_col5(:,1));
    final_il_range(:,2) = max(tg_avg_RMS_col5(:,1));
    
    starting_xl_range(:,1) = min(tg_avg_RMS_col4(:,2));
    starting_xl_range(:,2) = max(tg_avg_RMS_col4(:,2));
    final_xl_range(:,1) = min(tg_avg_RMS_col5(:,2));
    final_xl_range(:,2) = max(tg_avg_RMS_col5(:,2));   
end

% Looping through the dataset to select sub-sample (if any)
count = 1;
for ii = 1:length(tg_avg_RMS_col4)
    if tg_avg_RMS_col4(ii,1) >= starting_il_range(:,1) && tg_avg_RMS_col4(ii,1) <= starting_il_range(:,2) 
        if tg_avg_RMS_col4(ii,2) >= starting_xl_range(:,1) && tg_avg_RMS_col4(ii,2) <= starting_xl_range(:,2)
            tg_avg_starting_RMS(count,:) = tg_avg_RMS_col4(ii,:);
            count = count + 1;
        end
    end
end

count = 1;
for ii = 1:length(tg_avg_RMS_col5)
    if tg_avg_RMS_col5(ii,1) >= final_il_range(:,1) && tg_avg_RMS_col5(ii,1) <= final_il_range(:,2) 
        if tg_avg_RMS_col5(ii,2) >= final_xl_range(:,1) && tg_avg_RMS_col5(ii,2) <= final_xl_range(:,2)
            tg_avg_final_RMS(count,:) = tg_avg_RMS_col5(ii,:);
            count = count + 1;
        end
    end
end

% Calculating the histograms for the datsets and the statistics 
hist(tg_avg_RMS_col4(:,4),100)
xlabel('Window average trim sum value','FontSize', 15)
ylabel('Count','FontSize', 15)

hold all

hist(tg_avg_RMS_col5(:,4),100)
h = legend('Starting model','Final model');
title(sprintf('Average trim sum horizon histogram (il%d-%d xl%d-%d)',starting_il_range(:,1),starting_il_range(:,2),starting_xl_range(:,1), starting_xl_range(:,2), 'FontSize', 15);
% Setting the colours of the starting and final modal histograms to be
% different, and changing the legend accordingly
hline = findobj(gcf,'Type','patch');
set(hline(3), 'FaceColor', 'red', 'EdgeColor', 'red');
alpha(hline(3),.5)
xlim([0,10]);
set(h, 'FontSize', 15)
set(h, 'Color', 'none');
set(h, 'Box', 'off');

% Calculating the modal bin range for the starting model histogram
[n,xout] = hist(tg_avg_RMS_col4(:,4),100);
width = xout(2) - xout(1);
xout(1);
xout(2);
xoutmin = xout-width/2;
xoutmax = xout+width/2;

for jj = 1:length(n);
    if n(1,jj) == max(n)
      text(0.8,0.18, sprintf('Starting model modal bin range: %5.3f - %5.3f', xoutmin(1,jj), xoutmax(1,jj)),'fontsize',15,'Units', 'normalized')  
    end
end

% Calculating the modal bin range for the final model histogram
[n,xout] = hist(tg_avg_RMS_col5(:,4),100);
width = xout(2) - xout(1);
xout(1);
xout(2);
xoutmin = xout-width/2;
xoutmax = xout+width/2;

for kk = 1:length(n);
    if n(1,kk) == max(n)
      text(0.8,0.16, sprintf('Final model modal bin range: %5.3f - %5.3f', xoutmin(1,kk), xoutmax(1,kk)),'fontsize',15,'Units', 'normalized')  
    end
end

start_kurt = kurtosis(tg_avg_RMS_col4(:,4));
final_kurt = kurtosis(tg_avg_RMS_col5(:,4));
start_skew = skewness(tg_avg_RMS_col4(:,4));
final_skew = skewness(tg_avg_RMS_col5(:,4));

text(0.8,0.14, sprintf('Starting model kurtosis: %5.3f',start_kurt),'fontsize',15,'Units', 'normalized')
text(0.8,0.12, sprintf('Final model kurtosis: %5.3f',final_kurt),'fontsize',15,'Units', 'normalized')
text(0.8,0.10, sprintf('Starting model skewness: %5.3f',start_skew),'fontsize',15,'Units', 'normalized')
text(0.8,0.08, sprintf('Final model skewness: %5.3f',final_skew),'fontsize',15,'Units', 'normalized')
figure
% -------------------------------------------------------------------------
% min_il = min(tg_avg_RMS_col4(:,1));
% max_il = max(tg_avg_RMS_col4(:,1));
% 
% min_xl = min(tg_avg_RMS_col4(:,2));
% max_xl = max(tg_avg_RMS_col4(:,2));
% 
% pos_il = [min_il:2:max_il];
% pos_xl = [min_xl:1:max_xl];
% pos_il =  pos_il';
% pos_xl =  pos_xl';
% test_grid = meshgrid(pos_il,pos_xl);
% test_grid2 = meshgrid(pos_xl,pos_il);
% test_grid2 = test_grid2';
% gh = [test_grid(:) test_grid2(:)];
% [loca locb] = ismember(tg_avg_RMS_col4(:,1:2),gh,'rows');
% slice_col(locb) = tg_avg_RMS_col4(:,4);
% imagesc(reshape(slice_col,size(pos_xl,1),size(pos_il,1)))
% colorbar;
% title('Horizon slice of average summed trim values');
% axis off;
% figure
% 
% min_il = min(tg_avg_RMS_col5(:,1));
% max_il = max(tg_avg_RMS_col5(:,1));
% 
% min_xl = min(tg_avg_RMS_col5(:,2));
% max_xl = max(tg_avg_RMS_col5(:,2));
% 
% pos_il = [min_il:2:max_il];
% pos_xl = [min_xl:1:max_xl];
% pos_il =  pos_il';
% pos_xl =  pos_xl';
% test_grid = meshgrid(pos_il,pos_xl);
% test_grid2 = meshgrid(pos_xl,pos_il);
% test_grid2 = test_grid2';
% gh = [test_grid(:) test_grid2(:)];
% [loca locb] = ismember(tg_avg_RMS_col5(:,1:2),gh,'rows');
% slice_col(locb) = tg_avg_RMS_col5(:,4);
% imagesc(reshape(slice_col,size(pos_xl,1),size(pos_il,1)))
% colorbar;
% title('Horizon slice of average summed trim values');
% axis off;
% figure
end

%%
load window_length.mat;
no_files = size(window_length,1);

for ii = 1:no_files;
file_name = fullfile(sprintf('tg_avg_fault_%d_RMS_%d-%d', window_length(ii,3), window_length(ii,1), window_length(ii,2)));
% file_name = fullfile(sprintf('tg_avg_RMS'));
file = dlmread(file_name,'\t');

% Remove NANs
log2 = file(:,5) ~= 1.0000e+30;
log1 = file(:,4) ~= 1.0000e+30;

% Separating the starting and final model values (columns 4 & 5,
% respectively) 
tg_avg_RMS_col4 = file(log1,1:4);
tg_avg_RMS_col5 = [file(log2,1:3) file(log2,5)];

hist(tg_avg_RMS_col4(:,4),100)
xlabel('Window average trim sum value','FontSize', 15)
ylabel('Count','FontSize', 15)

hold all

hist(tg_avg_RMS_col5(:,4),100)
h = legend('Starting model','Final model');
a = title(sprintf('Fault %d average trim sum horizon histogram - window length %d-%dm',window_length(ii,3), window_length(ii,1), window_length(ii,2)));
set(a,'fontsize', 15);
hline = findobj(gcf,'Type','patch');
set(hline(3), 'FaceColor', 'red', 'EdgeColor', 'red');
% set(hline(4), 'EdgeColor', 'blue');
alpha(hline(3),.5)
xlim([0,10]);
set(h, 'FontSize', 15)
set(h, 'Color', 'none');
set(h, 'Box', 'off');

[n,xout] = hist(tg_avg_RMS_col4(:,4),100);
width = xout(2) - xout(1);
xout(1);
xout(2);
xoutmin = xout-width/2;
xoutmax = xout+width/2;

for jj = 1:length(n);
    if n(1,jj) == max(n)
      text(0.8,0.18, sprintf('Starting model modal bin range: %5.3f - %5.3f', xoutmin(1,jj), xoutmax(1,jj)),'fontsize',15,'Units', 'normalized')  
    end
end

[n,xout] = hist(tg_avg_RMS_col5(:,4),100);
width = xout(2) - xout(1);
xout(1);
xout(2);
xoutmin = xout-width/2;
xoutmax = xout+width/2;

for kk = 1:length(n);
    if n(1,kk) == max(n)
      text(0.8,0.16, sprintf('Final model modal bin range: %5.3f - %5.3f', xoutmin(1,kk), xoutmax(1,kk)),'fontsize',15,'Units', 'normalized')  
    end
end

start_kurt = kurtosis(tg_avg_RMS_col4(:,4));
final_kurt = kurtosis(tg_avg_RMS_col5(:,4));
start_skew = skewness(tg_avg_RMS_col4(:,4));
final_skew = skewness(tg_avg_RMS_col5(:,4));

text(0.8,0.14, sprintf('Starting model kurtosis: %5.3f',start_kurt),'fontsize',15,'Units', 'normalized')
text(0.8,0.12, sprintf('Final model kurtosis: %5.3f',final_kurt),'fontsize',15,'Units', 'normalized')
text(0.8,0.10, sprintf('Starting model skewness: %5.3f',start_skew),'fontsize',15,'Units', 'normalized')
text(0.8,0.08, sprintf('Final model skewness: %5.3f',final_skew),'fontsize',15,'Units', 'normalized')
figure
end

%%
load dome_window_length.mat;
no_files = size(dome_window_length,1);

for ii = 1:no_files;
file_name = fullfile(sprintf('tg_avg_dome_%d_RMS_%d-%d', dome_window_length(ii,3), dome_window_length(ii,1), dome_window_length(ii,2)));
% file_name = fullfile(sprintf('tg_avg_RMS'));
file = dlmread(file_name,'\t');

% Remove NANs
log2 = file(:,5) ~= 1.0000e+30;
log1 = file(:,4) ~= 1.0000e+30;

% Separating the starting and final model values (columns 4 & 5,
% respectively) 
tg_avg_RMS_col4 = file(log1,1:4);
tg_avg_RMS_col5 = [file(log2,1:3) file(log2,5)];

hist(tg_avg_RMS_col4(:,4),100)
xlabel('Window average trim sum value','FontSize', 15)
ylabel('Count','FontSize', 15)

hold all

hist(tg_avg_RMS_col5(:,4),100)
h = legend('Starting model','Final model');
a = title(sprintf('Dome %d average trim sum horizon histogram - window length %d-%dm',dome_window_length(ii,3), dome_window_length(ii,1), dome_window_length(ii,2)));
set(a,'fontsize', 15);
hline = findobj(gcf,'Type','patch');
set(hline(3), 'FaceColor', 'red', 'EdgeColor', 'red');
% set(hline(4), 'EdgeColor', 'blue');
alpha(hline(3),.5)
xlim([0,20]);
set(h, 'FontSize', 15)
set(h, 'Color', 'none');
set(h, 'Box', 'off');

[n,xout] = hist(tg_avg_RMS_col4(:,4),100);
width = xout(2) - xout(1);
xout(1);
xout(2);
xoutmin = xout-width/2;
xoutmax = xout+width/2;

for jj = 1:length(n);
    if n(1,jj) == max(n)
      text(0.8,0.18, sprintf('Starting model modal bin range: %5.3f - %5.3f', xoutmin(1,jj), xoutmax(1,jj)),'fontsize',15,'Units', 'normalized')  
    end
end

[n,xout] = hist(tg_avg_RMS_col5(:,4),100);
width = xout(2) - xout(1);
xout(1);
xout(2);
xoutmin = xout-width/2;
xoutmax = xout+width/2;

for kk = 1:length(n);
    if n(1,kk) == max(n)
      text(0.8,0.16, sprintf('Final model modal bin range: %5.3f - %5.3f', xoutmin(1,kk), xoutmax(1,kk)),'fontsize',15,'Units', 'normalized')  
    end
end

start_kurt = kurtosis(tg_avg_RMS_col4(:,4));
final_kurt = kurtosis(tg_avg_RMS_col5(:,4));
start_skew = skewness(tg_avg_RMS_col4(:,4));
final_skew = skewness(tg_avg_RMS_col5(:,4));

text(0.8,0.14, sprintf('Starting model kurtosis: %5.3f',start_kurt),'fontsize',15,'Units', 'normalized')
text(0.8,0.12, sprintf('Final model kurtosis: %5.3f',final_kurt),'fontsize',15,'Units', 'normalized')
text(0.8,0.10, sprintf('Starting model skewness: %5.3f',start_skew),'fontsize',15,'Units', 'normalized')
text(0.8,0.08, sprintf('Final model skewness: %5.3f',final_skew),'fontsize',15,'Units', 'normalized')
figure
end
close