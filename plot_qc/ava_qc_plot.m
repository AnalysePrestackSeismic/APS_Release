function [] =ava_qc_plot(job_meta_path,i_block)
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
% Plots AVA of RMS Amplitude,
% Plots AVA of Variance
%%
close all;
job_meta = load(job_meta_path);

plotype =2;% plot type 1 menas vertical plot 2 means horizontal plot

qc_filepath = strcat(job_meta.ava_qc_directory,'ava_qc_',i_block,'.mat');
ava_qc=load(qc_filepath);

% aa= isfield(job_meta,angle);    % check if the job meta file has reference to angles
% 
% if aa==1
    % usually for angle stack based scanning
%     angles=zeros(size(job_meta.angle,1),1);
%     for pp=1:size(job_meta.angle,1)
%         angles(pp)=mean(job_meta.angle{pp});
%     end
%     clear pp
% else
%     % usually for angle gather based scanning
    angles=job_meta.tkey_min:job_meta.tkey_inc:job_meta.tkey_max;
% end
% angles=sind(angles).^2;
figure(1);
subplot(1,2,1);plot(ava_qc.rms);
xlabel( 'angle'); ylabel( 'RMS Amplitude');
subplot(1,2,2);plot(ava_qc.var);
xlabel( 'angle'); ylabel( 'Variance');

% figure(2);
% for k=1:ava_qc.n_win
%     t=strcat('z below wb: ',num2str(((k-0.5)*job_meta.ns_overlap_qc*job_meta.s_rate/1000)),'m/ms');
%     subplot(ava_qc.n_win,2,(2*k-1));
%     scatter(angles,ava_qc.rms(:,k),'MarkerEdgeColor','b','MarkerFaceColor','c');xlim([min(angles) max(angles)]);
%     xlabel( 'angle'); ylabel( 'RMS Amp');grid on;
%     title (t);
% 
%     
%     subplot(ava_qc.n_win,2,2*k);
%     scatter(angles,ava_qc.var(:,k),'MarkerEdgeColor','r','MarkerFaceColor','y');xlim([min(angles) max(angles)]);
%     xlabel( 'angle'); ylabel( 'Variance');grid on;
%     title (t);
% end

figure(3);
for k=1:ava_qc.n_win
    t=strcat('z below wb: ',num2str(((k-0.5)*job_meta.ns_overlap_qc*job_meta.s_rate/1000)),'m/ms');
    if plotype==1
        subplot(ava_qc.n_win,2,(2*k-1));
    elseif plotype==2
        subplot(2,ava_qc.n_win,k);
    end
    scatter(angles,(ava_qc.rms(:,k)-min(ava_qc.rms(:,k)))/(max(ava_qc.rms(:,k))-min(ava_qc.rms(:,k))),'MarkerEdgeColor','b','MarkerFaceColor','c');
    xlim([min(angles) max(angles)]);
    xlabel( 'angle'); ylabel( 'Norm RMS');grid on;
    title (t);
    
    if plotype==1
        subplot(ava_qc.n_win,2,2*k);
    elseif plotype==2
        subplot(2,ava_qc.n_win,ava_qc.n_win+k);
    end
    
    scatter(angles,(ava_qc.var(:,k)-min(ava_qc.var(:,k)))/(max(ava_qc.var(:,k))-min(ava_qc.var(:,k))),'MarkerEdgeColor','r','MarkerFaceColor','y');xlim([min(angles) max(angles)]);
    xlabel( 'angle'); ylabel( 'Norm Var');grid on;
    title (t);
end


end