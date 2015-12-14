function [thickness_ap,top_amplitude]=tuning_curve(wavelet,tune_plot)
%% Function to transform a wavelet into a tuning curve

l=length(wavelet);
pad=zeros(1,l);
wavelet_padded=[pad wavelet pad];
ind1=l:(2*l-1);
wavelet2=-1*wavelet;
figure(5);
px=0;
r=floor(l/20);
top_index=zeros(1,r+1);
bot_index=zeros(1,r+1);
top_amplitude = zeros(1,r+1);
bot_amplitude = zeros(1,r+1);

for ind2=(-1*r):0
    x=wavelet2+wavelet_padded(ind1+ind2);
%     if tune_plot==1;
%         if px/20==floor(px/20) || ind2==0
%             plot(x);ylim([-2 2]);pause(0.1);
%         end
%     end
     px=px+1;
     %pick top and base
     [top_amplitude(px),top_index(px)]=max(x);
     [bot_amplitude(px),bot_index(px)]=min(x);
     
end

thickness_ap=  bot_index-top_index;
if tune_plot==1;
    figure(6);plot(thickness_ap,top_amplitude,'-r','LineWidth',3);
end
end
