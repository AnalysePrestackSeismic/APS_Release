[filter_files,nfiles] = directory_scan({'/data/TZA/dtect/2015_Kusini_pgs_trim/Misc/'},'crossplot.xyz');
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

[~,~,~,~,Z,~,~,~,ig,gradient_md,intercept_md,minenergy_md] = ... 
    importfile([filter_files.path{1} filter_files.names{6}]);

decim = floor(length(intercept_md)/n_points);
intercept_md = intercept_md(1:decim:end);
gradient_md = gradient_md(1:decim:end);
minenergy_md = minenergy_md(1:decim:end);

[polcoord(:,1),polcoord(:,2)] = cart2pol(intercept_md,gradient_md);
polcoord(:,3) = 1:1:length(intercept_md);
polcoord(:,1) = round(polcoord(:,1)*100)/100;
polcoord(:,2) = round(polcoord(:,2)*100)/100;

polcoord = sortrows(polcoord,[1 -2]);

thet = unique(polcoord(:,1));

count = 1;
for i_thet = 1:1:length(thet)    
    rad = polcoord(polcoord(:,1) == thet(i_thet),3);
    polout(count,:) = rad(1);
    count = count + 1;
end
smoothl = 3;
smoothcoef = ones(1, smoothl)/smoothl;
intercept_md = filter(smoothcoef,1,intercept_md(polout));
gradient_md = filter(smoothcoef,1,gradient_md(polout));

for i_file = 1:1:nfiles

[Inline1,Xline,~,~,Z,~,~,~,ig,gradient,intercept,minenergy] = ... 
    importfile([filter_files.path{i_file} filter_files.names{i_file}]);

% [polcoord(:,1),polcoord(:,2)] = cart2pol(intercept_md,gradient_md);
% polcoord(:,3) = 1:1:length(intercept_md);
% polcoord(:,1) = round(polcoord(:,1)*100)/100;
% polcoord(:,2) = round(polcoord(:,2)*100)/100;
% 
% polcoord = sortrows(polcoord,[1 -2]);
% 
% thet = unique(polcoord(:,1));
% 
% count = 1;
% for i_thet = 1:1:length(thet)    
%     rad = polcoord(polcoord(:,1) == thet(i_thet),3);
%     polout(count,:) = rad(1);
%     count = count + 1;
% end

% decimate
decim = floor(length(intercept)/n_points);
intercept = intercept(1:decim:end);
gradient = gradient(1:decim:end);
minenergy = minenergy(1:decim:end);
figure(1)
subplot(2,3,i_file); scatter(intercept,gradient,20,minenergy,'fill')
colormap(brewermap([],'*RdGy'))
caxis([-250 250])
colorbar
%subplot(2,3,i_file); scatter(intercept,gradient,20,'fill')
hold all
plot(intercept_md,gradient_md,'r','LineWidth',1.5)
hold off
legend(filter_files.names{i_file},'Location','NorthEast')
axis([-1500 1500 -1500 1500])
grid on

%hold all
%figure(1)
%scatter(intercept(1:decim:end,:),gradient(1:decim:end,:),5)

%figure(2)
%scatter3(intercept(1:decim:end,:),gradient(1:decim:end,:),minenergy(1:decim:end,:))

% [polcoord(:,1),polcoord(:,2)] = cart2pol(intercept,gradient);
% polcoord(:,3) = 1:1:length(intercept);
% polcoord(:,1) = round(polcoord(:,1)*100)/100;
% polcoord(:,2) = round(polcoord(:,2)*100)/100;
% 
% polcoord = sortrows(polcoord,[1 -2]);
% 
% thet = unique(polcoord(:,1));
% 
% count = 1;
% for i_thet = 1:1:length(thet)    
%     rad = polcoord(polcoord(:,1) == thet(i_thet),3);
%     polout(count,:) = rad(1);
%     count = count + 1;
% end
% 
% %scatter(intercept,gradient)
% hold all
% figure(1)
% scatter(intercept,gradient,10,'fill')
% 
% if i_file == nfiles
%     legend(filter_files.names,'Location','EastOutside')
%     axis([-2000 2000 -2000 2000])
% end
% 
% figure(2)
% patch(intercept(polout),gradient(polout),i_file,'FaceAlpha',.3)%,'EdgeColor',i_file)
% 
% if i_file == nfiles
%     legend(filter_files.names,'Location','EastOutside')
%     axis([-2000 2000 -2000 2000])
% end
% 
% figure(3)
% %scatter(intercept,gradient,10,'fill')
% %hold all
% plot(intercept(polout),gradient(polout))
% 
% if i_file == nfiles
%     legend(filter_files.names,'Location','EastOutside')
%     axis([-2000 2000 -2000 2000])
% end
% 
% figure(4)
% scatter(intercept,gradient,10,'fill')
% 
% if i_file == nfiles
%     legend(filter_files.names,'Location','EastOutside')
%     axis([-1000 0 -1000 0])
% end
% 
% hold off
%clear polcoord polout

end