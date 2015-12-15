function[]=ava_combinations(vp_shale,vs_shale,rho_shale,vp_sand,vs_sand,rho_sand,max_angle,freq_c,method)
%% Function to plot combinations of sands and shales as a matrix. Plots both the tuned ava matrix and the untured ava curve
close all;


n_shales=length(vp_shale);              % Number of shales
n_sands=length(vs_sand);                % Number of sands
angles=0:max_angle;                     % Array of Angle Values
n_comb =n_shales*n_sands;               % Number of Combinations
[temp1,temp2,temp3]=ava_wedge(vp_shale(1),vs_shale(1),rho_shale(1),vp_sand(1),vs_sand(1),rho_sand(1),max_angle,freq_c,method,0);
temp_int=temp1(:,1);
[~,index_tuning]= min(temp_int);        % Find index where maximum constructive interference occurs
size_ava_mat=size(temp1);               % calculate size of 2D ava matrices
ava_refl_comb=zeros(n_comb,size_ava_mat(1),size_ava_mat(2));    % pre allocate memory for storing the ava matrices for different combinations
thickness_true=zeros(n_comb,length(temp2));                     % pre allocate memory for storing the true thickness for different combinations
thickness_ap=zeros(n_comb,length(temp3));                       % pre allocate memory for storing the apparent thickness for different combinations
n_h=length(temp2);                                              % Number of horizontal samples( traces)
clear temp1 temp2 temp3 temp_int;                               % Clear redundant variables

% Calculate AVA for dirrent combinations of sands and shales
for j=1:n_sands
    for k=1:n_shales
        n=(j-1)*n_sands+k;
        [ava_refl_comb(n,:,:),thickness_true(n,:),thickness_ap(n,:)] =ava_wedge(vp_shale(k),vs_shale(k),rho_shale(k),vp_sand(j),vs_sand(j),rho_sand(j),max_angle,freq_c,method,0);
        
    end
end

r_max= max(max(max(ava_refl_comb)));    % Maximum reflectivity
r_min= min(min(min(ava_refl_comb)));    % Minimum Reflectivity
abs_max=max(max(max(abs(ava_refl_comb))));

%% Plot results

% Plot AVA matrix  with tuning
figure(11);
for j=1:n_sands
    for k=1:n_shales
        n=(j-1)*n_sands+k;
        subplot(n_sands,n_shales,n);
        time_thickness_true=thickness_true(n,:);
        imagesc([squeeze(ava_refl_comb(n,:,:))],[r_min r_max]);colorbar
        t=sprintf('Tuned AVA Shale: %d Sand: %d',k,j);
        title(t);
        xlabel('Angle (deg)');
        set(gca,'YTick',200:200:n_h);
        set(gca,'YTickLabel',time_thickness_true(200:200:n_h));
        xlabel('True Time Thickness (ms)');
        
    end
end

% Plot AVA curves when not tuned and when fully tuned
figure(12);
for j=1:n_sands
    for k=1:n_shales
        n=(j-1)*n_sands+k;
        subplot(n_sands,n_shales,n);
        plot(angles,squeeze(ava_refl_comb(n,end,:)),'LineWidth',3);
        hold on;
        plot(angles,squeeze(ava_refl_comb(n,index_tuning,:)),'--b','LineWidth',3);
        hold off;
        ylim([-abs_max abs_max]);xlim([0 max_angle]); 
        hold on;
        plot([0 max_angle],[0 0],'-g','LineWidth',2)
        t=sprintf('AVA Shale: %d Sand: %d',k,j);
        title(t);
        xlabel('Angle');
        ylabel('Reflectivity');
    end
end
   
        
end
