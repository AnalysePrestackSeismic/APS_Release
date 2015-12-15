function [] = wavelet_analysis(input_dir)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

cd(input_dir);
    
% Figure out the number of files in the current directory
[~,nfiles] = (system('ls -B | wc -l')); 
nfiles = str2double(nfiles);

% Read all filenames and convert from ascii to double
[~,fnames] = system('ls -B1');
numeric = double(fnames);

% Preallocate memory for some variables
count = 2;
fname_index = zeros(1,nfiles+1);

% Loop to separate out each file name from one long character string
for ij= 1:length(fnames)
    if numeric(1,ij) == 10
        fname_index(1,count) = ij;
        count = count+1;
    end
end

% Loop to read each file as ascii into a cell array
for ik=1:nfiles
    files_in{ik} = fnames(1,(fname_index(1,ik)+1):(fname_index(1,ik+1)-1));
end

fprintf('\nThe following files have been found and loaded:\n')
for il = 1:nfiles
    fprintf('File %d: %s\n',il,files_in{il});
    w(:,il) = dlmread(files_in{il});
end

[nt ~] = size(w);
t = nt*4;
pos_time = 0:4:(t-4)/2;
neg_time = -4:-4:-(t-4)/2;
time = [fliplr(neg_time) pos_time];
if nt ~= length(time)
    w = w(2:end,:);
    nt = nt-1;
end
freq = 0:125/floor(nt/2):125;

% Plots in time domain
figure(1)
set(1,'Units','inches','Position',[0 0 10 10]);
% plot(time,w,'LineWidth',3);
plot(time,w(:,1),'b','LineWidth',3);
legend(files_in,'interpreter','none')
axis tight
title('Time domain','FontSize',14)
ylabel('Amplitude','FontSize',14)
xlabel('Time (ms)','FontSize',14)
set(gca,'FontSize',14)
print('-r300','-dtiff','time_domain.tiff')

afw = abs(fft(w));
afw = afw(1:ceil((nt+1)/2),:);
ylim = -max(max(afw));
aoctaves = [125,125,125/2,125/2,125/4,125/4,125/8,125/8,125/16,125/16,125/32,125/32,125/64,125/64,125/128,125/128;5,ylim,5,ylim,5,ylim,5,ylim,5,ylim,5,ylim,5,ylim,5,ylim]';

% Plots in frequency domain
figure(2)
set(2,'Units','inches','Position',[0 0 10 10]);
plot(freq,bsxfun(@minus,afw,max(afw)),'LineWidth',3)
axis([0,125,ylim,5])
legend(files_in,'interpreter','none')
title('Frequency domain','FontSize',14)
ylabel('Amplitude (below maximum)','FontSize',14)
xlabel('Frequency (Hz)','FontSize',14)
hold all
for i=1:8
    plot(aoctaves(1+(2*(i-1)):2*i,1),aoctaves(1+(2*(i-1)):2*i,2),'black','linestyle','--');
end
set(gca,'FontSize',14)
print('-r300','-dtiff','frequency_domain_amplitude.tiff')

pfw = 20*log10(bsxfun(@ldivide,max(afw),afw));
poctaves = [125,125,125/2,125/2,125/4,125/4,125/8,125/8,125/16,125/16,125/32,125/32,125/64,125/64,125/128,125/128;5,-40,5,-40,5,-40,5,-40,5,-40,5,-40,5,-40,5,-40]';

figure(3)
set(3,'Units','inches','Position',[0 0 10 10]);
% plot(freq,pfw,'LineWidth',3)
plot(freq,pfw(:,1),'b','LineWidth',3)
axis([0,125,-40,5])
legend(files_in,'interpreter','none')
title('Frequency domain','FontSize',14)
ylabel('Power (dB below maximum amplitude)','FontSize',14)
xlabel('Frequency (Hz)','FontSize',14)
hold all
for i=1:8
    plot(poctaves(1+(2*(i-1)):2*i,1),poctaves(1+(2*(i-1)):2*i,2),'black','linestyle','--');
end
set(gca,'FontSize',14)
print('-r300','-dtiff','frequency_domain_power.tiff')

iafw = abs(fft(w,256));
iafw = iafw(1:129,:);
ipfw = 20*log10(bsxfun(@ldivide,max(iafw),iafw));
ifreq = 0:125/128:125;
ifreq= ifreq';

% Dominant Frequncy
df = ismember(ipfw,0);

% Bandwidth
b = ismember(round(ipfw/10)*10,0);

% Central Frequency
cf = cumsum(b);
cf = (bsxfun(@times,cf,b));
for i = 1:nfiles
    cf(:,i) = ismember(cf(:,i),round((max(cf(:,i))-1)/2));
end
cf = logical(cf);

% Bandwdith edges
eb = diff(b);
eb = [zeros(1,nfiles);eb];
eb = abs(eb);

% Bandwdith lower edge
eb1 = cumsum(eb);
eb1 = ismember(eb1,[1 2]);
eb1 = diff(eb1);
eb1 = [zeros(1,nfiles);eb1];
eb1 = logical(eb1);

% Bandwdith upper edge
eb2 = cumsum(eb);
eb2 = ismember(eb2,2);
eb2 = diff(eb2);
eb2 = [zeros(1,nfiles);eb2];
eb2 = logical(eb2);
ebdfcf = eb+df+cf;

larrow = '\leftarrow';
rarrow = '\rightarrow';

figure(4)
set(4,'Units','inches','Position',[0 0 20 20]);
clims = [-40 0];
imagesc(ifreq,(1:1:nfiles),ipfw',clims)
hold all
% h1 = imagesc(ifreq,(1:1:nfiles),-20*ones(size(ipfw')));
% set(h1,'AlphaData',df')
h2 = imagesc(ifreq,(1:1:nfiles),-30*ones(size(ipfw')));
set(h2,'AlphaData',cf')
h3 = imagesc(ifreq,(1:1:nfiles),-40*ones(size(ipfw')));
set(h3,'AlphaData',eb')
for i=1:nfiles
    % if ifreq(df(:,i)) >= ifreq(cf(:,i))
        % text(ifreq(df(:,i)),i,sprintf('%s%.1f',larrow,ifreq(df(:,i))),'Interpreter','Tex','VerticalAlignment','Middle','HorizontalAlignment','Left','FontSize',14,'Color','white')
        % text(ifreq(cf(:,i)),i,sprintf('%.1f%s',ifreq(cf(:,i)),rarrow),'Interpreter','Tex','VerticalAlignment','Middle','HorizontalAlignment','Right','FontSize',14,'Color','white')
    % else
        % text(ifreq(df(:,i)),i,sprintf('%.1f%s',ifreq(df(:,i)),rarrow),'Interpreter','Tex','VerticalAlignment','Middle','HorizontalAlignment','Right','FontSize',14,'Color','white')
        text(ifreq(cf(:,i)),i,sprintf('%s%.1f',larrow,ifreq(cf(:,i))),'Interpreter','Tex','VerticalAlignment','Middle','HorizontalAlignment','Left','FontSize',14,'Color','white')
    % end
    text(ifreq(eb1(:,i)),i,sprintf('%s%.1f',larrow,ifreq(eb1(:,i))),'VerticalAlignment','Middle','HorizontalAlignment','Left','FontSize',14,'Color','white')
    text(ifreq(eb2(:,i)),i,sprintf('%s%.1f',larrow,ifreq(eb2(:,i))),'VerticalAlignment','Middle','HorizontalAlignment','Left','FontSize',14,'Color','white')
end
title('Spectrum Matrix','FontSize',14)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Wavelet','FontSize',14)
ylabel(colorbar,'Power (dB below maximum amplitude)','FontSize',14)
set(gca,'YTick',(1:length(files_in)))
set(gca,'yTickLabel',files_in)
set(gca,'FontSize',14)
print('-r300','-dtiff','spectrum_matrix.tiff')

end

