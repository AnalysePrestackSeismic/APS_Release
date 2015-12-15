clear;

file1 = '1_sliceopt_cgsh_grad_in300_xline650_wheeler' 
file2 = '1_sliceopt_cgsh_int_in300_xline650_wheeler'
file3 = '101_sliceopt_cgsh_data_EEI_linear_chi'
timefile = '1_sliceopt_cgsh_time_intgrad_in300_xline650_wheeler'
nt = 466; 
[backgroundtrend] = crossplotter6(file1,file2,466,1,4,466,220,650,300,1,1)
file4 = '220_EEI_data'
fid1 = fopen('1_sliceopt_cgsh_int_in300_xline650_wheeler');
fid2 = fopen('1_sliceopt_cgsh_grad_in300_xline650_wheeler');
fid3 = fopen('101_sliceopt_cgsh_data_EEI_linear_chi');
fid4 = fopen('220_EEI_data');

% 
% for g = 1:466
%     data1(:,1) = fread(fid1,195000,'float32');
%     data2(:,1) = fread(fid2,195000,'float32');
%     data3(:,1) = fread(fid3,195000,'float32');
%     data4(:,1) = fread(fid4,195000,'float32');
%     
%     for j = 1:length(data1)
%         
%         if data1(j,1) > 1e29;
%             data1(j,1) = NaN;
%         end
%            if data2(j,1) > 1e29;
%             data2(j,1) = NaN;
%         end     
%           if data3(j,1) > 1e29;
%             data3(j,1) = NaN;
%         end      
%            if data4(j,1) > 1e29;
%             data4(j,1) = NaN;
%         end     
%         
%     end
%     
%             
%     [prob1,pos1] = ecdf(data1);
%     [prob2,pos2] = ecdf(data2);
%     [prob3,pos3] = ecdf(data3);
%     [prob4,pos4] = ecdf(data4);
%             
%     for k = 1:99
%     b = k*0.01;
%     n5(g,k) = findnearest(b,prob1);
%     n6(g,k) = findnearest(b,prob2);
%     n7(g,k) = findnearest(b,prob3);
%     n8(g,k) = findnearest(b,prob4);
%     
%     val5(g,k) = pos1(n5(g,k),1);
%     val6(g,k) = pos2(n6(g,k),1);
%     val7(g,k) = pos3(n7(g,k),1);
%     val8(g,k) = pos4(n8(g,k),1);
%     
%     
%     end
%     
% end


[average1,stddev1, val1] = sas_stats_edit(file1,timefile,4,466,50,0.5,201)
[average2,stddev2, val2] = sas_stats_edit(file2,timefile,4,466,50,0.5,202)
[average3,stddev3, val3] = sas_stats_edit(file3,timefile,4,466,50,0.5,204)
[average4,stddev4, val4] = sas_stats_edit(file4,timefile,4,466,50,0.5,205)

x = -60:0.2:60;

for j = 1:466

max1 = (1/(stddev1(j)*sqrt(2*pi)));
max2 = (1/(stddev2(j)*sqrt(2*pi)));
max3 = (1/(stddev3(j)*sqrt(2*pi)));
max4 = (1/(stddev4(j)*sqrt(2*pi)));   

norm_gaus1(j,:) = ((1/(stddev1(j)*sqrt(2*pi)))*exp(-0.5*(((x-average1(j))/(stddev1(j))).^2)))/max1;
norm_gaus2(j,:) = ((1/(stddev2(j)*sqrt(2*pi)))*exp(-0.5*(((x-average2(j))/(stddev2(j))).^2)))/max2;
norm_gaus3(j,:) = ((1/(stddev3(j)*sqrt(2*pi)))*exp(-0.5*(((x-average3(j))/(stddev3(j))).^2)))/max3;
norm_gaus4(j,:) = ((1/(stddev4(j)*sqrt(2*pi)))*exp(-0.5*(((x-average4(j))/(stddev4(j))).^2)))/max4;

end


% x_axis = 0.01:0.01:0.99;
% 
% 
% figure('Color',[1 1 1])
% subplot(2,1,1)
% imagesc(val2,'XData',x_axis,[-60 60])
% c = colorbar
% set(c,'FontSize',18)
% ylabel(c,'Intercept data')
% title({'z - dependant Intercept probability density function using cdf in wheeler domain'},...
%     'FontSize',18)
% ylabel({'Slice number'},'FontSize',20);
% set(gca,'XTick',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9])
% subplot(2,1,2)
% imagesc(val5,'XData',x_axis,[-60 60])
% c = colorbar
% set(c,'FontSize',18)
% ylabel(c,'Intercept data value')
% title({'z - dependant Intercept probability density function  using ecdf in wheeler domain'},...
%    'FontSize',18)
% ylabel({'Slice number'},'FontSize',20);
% set(gca,'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
% 

% figure('Color',[1 1 1])
% subplot(2,1,1)
% imagesc(val1,'XData',x_axis,[-80 80])
% c = colorbar
% set(c,'FontSize',18)
% ylabel(c,'Gradient data')
% title({'z - dependant Gradient probability density function using cdf in wheeler domain'},...
%     'FontSize',18)
% ylabel({'Slice number'},'FontSize',20);
% set(gca,'XTick',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9])
% subplot(2,1,2)
% imagesc(val6,'XData',x_axis,[-80 80])
% c = colorbar
% set(c,'FontSize',18)
% ylabel(c,'Gradient data value')
% title({'z - dependant Gradient probability density function  using ecdf in wheeler domain'},...
%    'FontSize',18)
% ylabel({'Slice number'},'FontSize',20);
% set(gca,'XTick',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9])
% %%
% figure('Color',[1 1 1])
% 
% subplot(2,1,1)
% imagesc(val3,'XData',x_axis,[-80 80])
% c = colorbar
% set(c,'FontSize',18)
% ylabel(c,'EEI data- Linear background trend')
% title({'z - dependant EEI- Linear background trend probability density function using cdf in wheeler domain'},...
%     'FontSize',18)
% ylabel({'Slice number'},'FontSize',20);
% set(gca,'XTick',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9])
% 
% subplot(2,1,2)
% imagesc(val7,'XData',x_axis,[-80 80])
% c = colorbar
% set(c,'FontSize',18)
% ylabel(c,'EEI Data - Linear background trend')
% title({'z - dependant EEI- Linear background trend probability density function  using ecdf in wheeler domain'},...
%    'FontSize',18)
% ylabel({'Slice number'},'FontSize',20);
% set(gca,'XTick',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9])
% 
% 
% %%
% 
% 
% 
% figure('Color',[1 1 1])
% subplot(2,1,1)
% imagesc(val4,'XData',x_axis,[-80 80])
% c = colorbar
% set(c,'FontSize',18)
% ylabel(c,'EEI data - Principle component analysis background trend')
% title({'z - dependant EEI - Principle component analysis background trend probability density function using cdf in wheeler domain'},...
%     'FontSize',18)
% ylabel({'Slice number'},'FontSize',20);
% set(gca,'XTick',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9])
% subplot(2,1,2)
% imagesc(val8,'XData',x_axis,[-80 80])
% c = colorbar
% set(c,'FontSize',18)
% ylabel(c,'EEI Data- Principle component analysis background trend')
% title({'z - dependant EEI- Principle component analysis background trend probability density function  using ecdf in wheeler domain'},...
%    'FontSize',18)
% ylabel({'Slice number'},'FontSize',20);
% set(gca,'XTick',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9])
% 
figure('Color',[1 1 1])
imagesc(norm_gaus2,'XData',x);
title({'z-dependant gaussian models for wheeler flattened'...
    'intercept'},...
    'FontSize',30);
xlabel({'Intercept model Value'},'FontSize',30);
ylabel({'Slice number'},'FontSize',30);
set(gca,'XTick',[-60 -40 -20 0 20 40 60],'FontSize',30)
c = colorbar
set(c,'FontSize',20)
ylabel(c,'Normalised Probability Density')

figure('Color',[1 1 1])
imagesc(norm_gaus1,'XData',x);
title({'z-dependant gaussian models for wheeler flattened'...
    'gradient'},...
    'FontSize',30);
xlabel({'Gradient model value'},'FontSize',30);
ylabel({'Slice number'},'FontSize',30);
set(gca,'XTick',[-60 -40 -20 0 20 40 60],'FontSize',30)
c = colorbar
set(c,'FontSize',20)
ylabel(c,'Normalised Probability Density')

figure('Color',[1 1 1])
imagesc(norm_gaus4,'XData',x);
title({'z-dependant gaussian models for wheeler flattened EEI volume, generated using slice dependant PCA \chi_{0}'},...
    'FontSize',30);
xlabel({'EEI model value'},'FontSize',30);
ylabel({'Slice number'},'FontSize',30);
set(gca,'XTick',[-60 -40 -20 0 20 40 60],'FontSize',30)
c = colorbar
set(c,'FontSize',20)
ylabel(c,'Normalised Probability Density')

figure('Color',[1 1 1])
imagesc(norm_gaus3,'XData',x);
title({'z-dependant gaussian models for wheeler flattened EEI volume, generated using linear \chi_{0}'},...
    'FontSize',30);
xlabel({'EEI model value'},'FontSize',30);
ylabel({'Slice number'},'FontSize',30);
set(gca,'XTick',[-60 -40 -20 0 20 40 60],'FontSize',30)
c = colorbar
set(c,'FontSize',20)
ylabel(c,'Normalised Probability Density')
% 
figure('Color',[1 1 1])
subplot(2,2,1)
imagesc(norm_gaus2,'XData',x);
title({'z-dependant gaussian models for wheeler flattened',...
    'intercept'},...
    'FontSize',30);
xlabel({'Intercept model Value'},'FontSize',30);
ylabel({'Slice number'},'FontSize',30);
set(gca,'XTick',[-60 -40 -20 0 20 40 60],'FontSize',30)
c = colorbar
set(c,'FontSize',20)
ylabel(c,'Normalised Probability Density')
%
subplot(2,2,2)
imagesc(norm_gaus1,'XData',x);
title({'z-dependant gaussian models for wheeler flattened'...
    'gradient'},...
    'FontSize',30);
xlabel({'Gradient model value'},'FontSize',30);
ylabel({'Slice number'},'FontSize',30);
set(gca,'XTick',[-60 -40 -20 0 20 40 60],'FontSize',30)
c = colorbar
set(c,'FontSize',20)
ylabel(c,'Normalised Probability Density')

subplot(2,2,3)
imagesc(norm_gaus4,'XData',x);
title({'z-dependant gaussian models for wheeler flattened'...
    'EEI volume generated using slice dependant PCA \chi_{0}'},...
    'FontSize',30);
xlabel({'EEI model value'},'FontSize',30);
ylabel({'Slice number'},'FontSize',30);
set(gca,'XTick',[-60 -40 -20 0 20 40 60],'FontSize',30)
c = colorbar
set(c,'FontSize',20)
ylabel(c,'Normalised Probability Density')
subplot(2,2,4)
imagesc(norm_gaus3,'XData',x);
title({'z-dependant gaussian models for wheeler flattened'...
    'EEI volume generated using linear \chi_{0}'},...
    'FontSize',30);
xlabel({'EEI model value'},'FontSize',30);
ylabel({'Slice number'},'FontSize',30);
set(gca,'XTick',[-60 -40 -20 0 20 40 60],'FontSize',30)
c = colorbar
set(c,'FontSize',20)
ylabel(c,'Normalised Probability Density')
%





