function[]=ava_wedge_depth_trends(trend_file_overburden,trend_file_reservoir,job_meta_path,wavelet_file,method)
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
% inputs: 
%     trend_file_overburden : path of trend las file for overburden
%     trend_file_reservoir : path of trend las file for reservoir
%     job_meta_path : path of job_meta file of digi_project
%     wavelet_file : name of wavlet file from digi project
%     method: 1 for avo 3 term aki richard's , 2 for zeopritz
% ouput:
%     void
% Displays:
%   the Vp, Vs and density trends
%   AVO based plots for specified depths
% Note:
%   The trend files should be las files with columns of depth Vp Vs rho in  that order
%%
depths=1000:1000:4000;
n_p=length(depths);
close all;
trends_overburden = read_las_file(trend_file_overburden);
trends_reservoir=read_las_file(trend_file_reservoir);
figure(1);
for kk=1:3
    subplot(1,3,kk);
    plot(trends_overburden.curves(:,kk+1),trends_overburden.curves(:,1),'-r','LineWidth',2);hold on;
    plot(trends_reservoir.curves(:,kk+1),trends_reservoir.curves(:,1),'-b','LineWidth',2);
    
    ylabel('depth')
    set(gca,'Ydir','reverse');
    switch kk
        case 1
        xlabel('vp'); title('Compressional Vel')
        plot([1200 5000],[depths;depths],'--g','LineWidth',3);
        xlim([1200 5000]);
        case 2
        xlabel('vs'); title('Shear Vel')
        plot([0 3000],[depths;depths],'--g','LineWidth',3);
        xlim([0 3000]);
        case 3
        xlabel('rhob'); title('Density')
        plot([1 3.3],[depths;depths],'--g','LineWidth',3);
        xlim([1 3.3]);
    end
end
hold off;


% for a=1:n_p
%   output{a}=ava_wedge_fd_digi(trend_file_overburden,trend_file_reservoir,job_meta_path,wavelet_file,depths(a),method,0);  
% end

figure(2);
for ii=1:n_p
    n_c=6;
    k=n_c*(ii-1);
    output_temp=ava_wedge_fd_digi(trend_file_overburden,trend_file_reservoir,job_meta_path,wavelet_file,depths(ii),method,0);  
    j=1;
    subplot(n_p,n_c,k+j);
    plot(output_temp.ang,output_temp.r_ava,'-g','LineWidth',3);hold on;
    plot([0 max(output_temp.ang)],[0 0],'k','LineWidth',2);hold off;
    xlabel('Angle');ylabel('Reflectivity');title('ava plot');
    ylim([-1 1]);
    j=j+1;
    
    subplot(n_p,n_c,k+j);
    plot(output_temp.t_wavelet,output_temp.wavelet_set,'LineWidth',2);
    xlabel('time');title('Wavlelet');
    xlim([0 max(output_temp.t_wavelet)]);
    j=j+1;
    
    subplot(n_p,n_c,k+j);
    imagesc(output_temp.ava_refl_top);
    n=size(output_temp.ava_refl_top,1);
    set(gca,'YTick',floor(n/10):floor(n/10):n);
    set(gca,'YTickLabel',output_temp.time_thickness_true(floor(n/10):floor(n/10):n));
    set(gca,'XTick',1:size(output_temp.ava_refl_top,2));
    set(gca,'XTickLabel',output_temp.ang(1:size(output_temp.ava_refl_top,2)));
    ylabel('True Time Thickness (ms)');   
    title('AVA Top');xlabel('Angle');
    j=j+1;
    
    subplot(n_p,n_c,k+j);
    plot(output_temp.ang,output_temp.ava_refl_top);hold on; 
    xlabel('Angle');ylabel('Reflectivity');title('ava plot reconstructed');
    plot([0 max(output_temp.ang)],[0 0],'k','LineWidth',2);
    ylim([-1 1]);hold off;
    j=j+1;
    
    subplot(n_p,n_c,k+j);
    scatter(output_temp.int,output_temp.grad,5,[0 1 1]);
    xlim([-max(abs(1.1*output_temp.int)) max(1.1*abs(output_temp.int))]);
    ylim([-max(abs(1.1*output_temp.grad)) max(1.1*abs(output_temp.grad))]);
    grid on;
    xlabel('intercept');ylabel('Gradient');
    title('Tuning IG Hodogram');
    j=j+1;
    
    subplot(n_p,n_c,k+j);
    scatter(output_temp.time_thickness_true,output_temp.chi_angle,5,[0 0 1]);
    xlim([0 max(abs(output_temp.time_thickness_true))]);
    ylim([0 360]); xlabel('True time thickness');
    ylabel('Polar_Angle(I-G)');
    
    
end

end