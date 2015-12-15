% program to plot EER distributions

function [] = histogram_plots()
%
%target area
first_samp = 400;           % 490;
last_samp = 700;           % 570;
first_trace = 200;           % 600;
last_trace = 1800;           % 1200;
%first_samp = 400;
%last_samp = 660;
%first_trace = 400;
%last_trace = 2200;
minhist = -1000;
maxhist = 1000;
minyhist = 0;
maxyhist = 4250;
maxystdhist = 180;
histcentres = minhist:10:maxhist;


off_filename='/data/TZA/segy/2013_pweza_chewa_presdm_pgs/gather_for_eage/EER_dtect_export/eer_digi_w_range_7_5_35_off_trim.sgy'; 
noq_filename='/data/TZA/segy/2013_pweza_chewa_presdm_pgs/gather_for_eage/EER_dtect_export/eer_digi_w_range_8_1_30_noQ_trim.sgy';
q_filename='/data/TZA/segy/2013_pweza_chewa_presdm_pgs/gather_for_eage/EER_dtect_export/eer_digi_w_range_7_5_35_Q_trim.sgy';

[off_meta off_ilxl_bytes off_traces]=segy_to_mat('189','193',off_filename);                                 % Load the intercept file for the block
[noq_meta noq_ilxl_bytes noq_traces]=segy_to_mat('189','193',noq_filename);                             % Load the gradient file for the block
[q_meta q_ilxl_bytes q_traces]=segy_to_mat('189','193',q_filename); 

%%

figure(1); imagesc(noq_traces); colormap(gray); caxis([-400 400]); title('EER section from no q angle gathers');                % show example data section

sub_off = off_traces(first_samp:last_samp,first_trace:last_trace);
sub_noq= noq_traces(first_samp:last_samp,first_trace:last_trace);
sub_q = q_traces(first_samp:last_samp,first_trace:last_trace);

figure(2); imagesc(sub_noq); colormap(gray); caxis([-400 400]); title('zoom area of EER from no q angle gathers');              % zoom into area of interest

%%

sub_off=sub_off(:);
sub_noq=sub_noq(:);
sub_q=sub_q(:);

offstd=std(sub_off);
noqstd=std(sub_noq);
qstd=std(sub_q);

%%

figure(3);              % Histograms for each dataset plotted separately
subplot(3,1,1); 
hist(sub_off,histcentres);
xlim([minhist maxhist]);
ylim([minyhist maxyhist]);
title('histogram for offset EER');
subplot(3,1,2); 
hist(sub_noq,histcentres);
xlim([minhist maxhist]);
ylim([minyhist maxyhist]);
title('histogram for no Q EER');
subplot(3,1,3); 
hist(sub_q,histcentres);
xlim([minhist maxhist]);
ylim([minyhist maxyhist]);
title('histogram for Q EER');


figure(4);              % Histograms for each dataset plotted separately with data within 2 standard deviations from the mean taken out
subplot(3,1,1); 
hist(sub_off(or(sub_off>(offstd*2),sub_off<(offstd*-2))),histcentres);
xlim([minhist maxhist]);
ylim([minyhist maxystdhist]);
title(strcat('histograms > or < 2*SD (',num2str(offstd),') for offset EER, norm. kurtosis =',num2str(kurtosis(sub_off)-3),' skewness =',num2str(skewness(sub_off)) ));
subplot(3,1,2); 
hist(sub_noq(or(sub_noq>(noqstd*2),sub_noq<(noqstd*-2))),histcentres);
xlim([minhist maxhist]);
ylim([minyhist maxystdhist]); 
title(strcat('histograms > or < 2*SD (',num2str(noqstd),') for no Q EER, norm. kurtosis =',num2str(kurtosis(sub_noq)-3),' skewness =',num2str(skewness(sub_noq)) ));
%title(strcat('histograms > or < 2*SD (',num2str(noqstd),') for no Q EER'));
subplot(3,1,3); 
hist(sub_q(or(sub_q>(qstd*2),sub_q<(qstd*-2))),histcentres);
xlim([minhist maxhist]);
ylim([minyhist maxystdhist]);
title(strcat('histograms > or < 2*SD (',num2str(qstd),') for Q EER, norm. kurtosis =',num2str(kurtosis(sub_q)-3),' skewness =',num2str(skewness(sub_q)) ));
%title(strcat('histograms > or < 2*SD (',num2str(qstd),') for Q EER'));


[freq_off,out_off]=hist(sub_off,histcentres);
[freq_noq,out_noq]=hist(sub_noq,histcentres);
[freq_q,out_q]=hist(sub_q,histcentres);
% figure(5);                                      % Histograms for each dataset plotted separately
% subplot(3,1,1);
% bar(out_off,freq_off./max(freq_off));
% ylim([0 1]);
% title('normalised EER distribution off');
% subplot(3,1,2);
% bar(out_noq,freq_noq./max(freq_noq));
% ylim([0 1]);
% title('normalised EER distribution noq');
% subplot(3,1,3);
% bar(out_q,freq_q./max(freq_q));
% ylim([0 1]);
% title('normalised EER distribution q');

%[freq_off,out_off]=hist(sub_off,histcentres);
%[freq_noq,out_noq]=hist(sub_noq,histcentres);
%[freq_q,out_q]=hist(sub_q,histcentres);
% figure(5);
% subplot(3,1,1);
% bar(out_off,freq_off/sum(freq_off)*10);
% ylim([0 1]);
% title('normalised EER distribution off');
% subplot(3,1,2);
% bar(out_noq,freq_noq/sum(freq_noq)*10);
% ylim([0 1]);
% title('normalised EER distribution noq');
% subplot(3,1,3);
% bar(out_q,freq_q/sum(freq_q)*10);
% ylim([0 1]);
% title('normalised EER distribution q');

% Normalise by peak
figure (5)
plot(out_off,freq_off./max(freq_off));
ylim([0 1]);
title('normalised EER distribution off');
hold all
plot(out_noq,freq_noq./max(freq_noq));
ylim([0 1]);
title('normalised EER distribution noq');
plot(out_q,freq_q./max(freq_q));
ylim([0 1]);
title('normalised EER distribution q');
hold off
legend('Angle Stacks','Angle Gathers','Angle-Q Gathers')
title('Distributions normalised by peak value')

% [freq_off,out_off]=hist(sub_off(or(sub_off>(offstd*2),sub_off<(offstd*-2))),histcentres);
% [freq_noq,out_noq]=hist(sub_noq(or(sub_noq>(noqstd*2),sub_noq<(noqstd*-2))),histcentres);
% [freq_q,out_q]=hist(sub_q(or(sub_q>(qstd*2),sub_q<(qstd*-2))),histcentres);
% figure(6);
% subplot(3,1,1);
% bar(out_off,freq_off./max(freq_off));
% ylim([0 1]);
% title('normalised EER distribution off');
% subplot(3,1,2);
% bar(out_noq,freq_noq./max(freq_noq));
% ylim([0 1]);
% title('normalised EER distribution noq');
% subplot(3,1,3);
% bar(out_q,freq_q./max(freq_q));
% ylim([0 1]);
% title('normalised EER distribution q');

% Normalise by area

[sub_off_norm,no_off] = histnorm(sub_off,100);
[sub_noq_norm,no_noq] = histnorm(sub_noq,100);
[sub_q_norm,no_q] = histnorm(sub_q,100);
figure (6)
plot(no_off,sub_off_norm)
hold all
plot(no_noq,sub_noq_norm)
plot(no_q,sub_q_norm)
legend('Angle Stacks','Angle Gathers','Angle-Q Gathers')
title('Distributions normalised by area')
hold off
% peak normalise


% figure (7)
% data = [sub_off,sub_noq,sub_q];
% hist(data,100)
% 
% figure
% h1 = hist(sub_off,histcentres);
% hold all
% h2 = hist(sub_noq,histcentres);
% 
% h3 = hist(sub_q,histcentres);
% hold off
% [freq_off,out_off]=hist(sub_off(or(sub_off>(offstd*2),sub_off<(offstd*-2))),histcentres);
% [freq_noq,out_noq]=hist(sub_noq(or(sub_noq>(noqstd*2),sub_noq<(noqstd*-2))),histcentres);
% [freq_q,out_q]=hist(sub_q(or(sub_q>(qstd*2),sub_q<(qstd*-2))),histcentres);
% figure(6);
% subplot(3,1,1);
% bar(out_off,freq_off/sum(freq_off)*10);
% ylim([0 1]);
% title('normalised EER distribution off');
% subplot(3,1,2);
% bar(out_noq,freq_noq/sum(freq_noq)*10);
% ylim([0 1]);
% title('normalised EER distribution noq');
% subplot(3,1,3);
% bar(out_q,freq_q/sum(freq_q)*10);
% ylim([0 1]);
% title('normalised EER distribution q');



% figure(5);
% subplot(3,1,1);
% kurtosis(sub_off)
% %moment(sub_off,4);
% %skewness(hist(sub_off,centres));
% subplot(3,1,2)
% kurtosis(sub_noq,histcentres);
% subplot(3,1,3);
% kurtosis(sub_q,histcentres);

end





