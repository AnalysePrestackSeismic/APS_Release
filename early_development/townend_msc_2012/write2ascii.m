fid='output/gullris_hough.txt';
cube_out=Lcube_hough(:,:,:);
start_time=1682;%0;%1682;   %ms
sample_rate=2;%1;%2;  %ms
nsamples=size(cube_out,1);
start_iline=2490;%1;%2490;
start_xline=1804;%1;%1804;

% write header
dlmwrite(fid,[start_time,sample_rate,nsamples],'delimiter',' ','newline','unix')

% write every trace
for i=1:size(cube_out,2)
    [size(cube_out,2),i]
    tic;
    for j=1:size(cube_out,3)
        dlmwrite(fid,[start_iline+i-1,start_xline+j-1,cube_out(:,i,j)'],...
            '-append','delimiter',' ','newline','unix');
    end
    if i==1,tot_time=toc*size(cube_out,2)/60;end
    disp(tot_time)
    tot_time=tot_time-toc/60;
end