function [output] = minivolread(fhandle,filenames,bitsperbyte,nil,nxl,nt,nil_stepout,nxl_stepout,nt_stepout,il_olap,xl_olap,t_olap)

% Bits to bytes conversion factor
%bitsperbyte = 4;

% Numbers of il/xl/t in whole volume
% nil = 301;
% nxl = 651;
% nt = 838;

% Sample increments of il/xl/t (unused right now, but needed if user wants to enter il/xl/t ranges)
il_inc = 1;
xl_inc = 2;
t_inc = 4;

% Stepouts for computation
% nil_stepout = 25;
% nxl_stepout = 25;
% nt_stepout = 50;

% Fractional overlap of computation blocks
% il_olap = 0.50;
% xl_olap = 0.50;
% t_olap = 0.50;

[nfiles, ~] = size(filenames);

% Number of il/xl/t to read per computation block
nil_read = (2*nil_stepout)+1;
nxl_read = (2*nxl_stepout)+1;
nt_read = (2*nt_stepout)+1;

% Number of il/xl/t to shift computation block
nil_shift = round(il_olap*nil_read);
nxl_shift = round(xl_olap*nxl_read);
nt_shift = round(t_olap*nt_read);

% Total number of times the block moves in il/xl/t directions
nil_block = floor(nil/(nil_read*il_olap));
nxl_block = floor(nxl/(nxl_read*xl_olap));
nt_block = floor(nt/(nt_read*t_olap));

% Number of il/xl/t to ignore at volume edges to avoid edge effects
nil_crop = nil-nil_block*nil_shift;
nxl_crop = nxl-nxl_block*nxl_shift;
nt_crop = nt-nt_block*nt_shift;

% Some loops to ensure volume to process can be centred within total volume
while 2*floor(nil_crop/2) ~= nil_crop
    nil_shift = nil_shift-1;
    nil_crop = nil-nil_block*nil_shift;
    if nil_shift < 1
        error('Cannot crop input volume evenly with the choosen stepout and overlap to avoid edge effects in the inline direction. Please alter inline stepout and/or overlap.')
    end
end

while 2*floor(nxl_crop/2) ~= nxl_crop
    nxl_shift = nxl_shift-1;
    nxl_crop = nxl-nxl_block*nxl_shift;
    if nxl_shift < 1
        error('Cannot crop input volume evenly with the choosen stepout and overlap to avoid edge effects in the crossline direction. Please alter crossline stepout and/or overlap.')
    end
end

while 2*floor(nt_crop/2) ~= nt_crop
    nt_shift = nt_shift-1;
    nt_crop = nt-nt_block*nt_shift;
    if nt_shift < 1
        error('Cannot crop input volume evenly with the choosen stepout and overlap to avoid edge effects in the Z direction. Please alter Z stepout and/or overlap.')
    end
end

% Work out starting il/xl/t
nil_start = (nil_crop/2)+1;
nxl_start = (nxl_crop/2)+1;
nt_start = (nt_crop/2)+1;

% Work out ending il/xl/t
nil_end = nil-((nil_crop/2)+1);
nxl_end = nxl-((nxl_crop/2)+1);
nt_end = nt-((nt_crop/2)+1);

% Point to the input file
for l = 1:nfiles
    fid{l} = fopen(filenames{l});
end

% Preallocate some variables
data{nfiles} = [];
output = zeros(nt_block-1,nxl_block-1,nil_block-1);

spmd(2)
% Data read and computation loopy loops
for p = 1:nt_block-1 % Move block in t direction
    for q = 1:nxl_block-1 % Move block in xl direction
        for r = 1:nil_block-1 % Move block in il direction
            % Alter byte starting location (if you are looking at an inline with xline increasing left to right, this is the top front left corner of the block)
            bytestart = bitsperbyte*(((nil_start+((r-1)*nil_shift)-1)*nxl*nt)+((nxl_start+((q-1)*nxl_shift)-1)*nt)+(nt_start+((p-1)*nt_shift)-1));
            % Move to start location
            fseek(fid{1},bytestart,'bof'); 
            % Loops to read the data
            for j = 1:nil_read
                for k = 1:nxl_read
                    for l = 1:nfiles
                        data{l}(:,k,j) = fread(fid{l},nt_read,'float32'); % Read some data at last
                    end
                    fseek(fid{1},bitsperbyte*(nt-nt_read),'cof'); % Skip to next t sample to read
                end
                fseek(fid{1},bitsperbyte*nt*(nxl-nxl_read),'cof'); % Skip to next inline to read traces from
            end
            [nt_data,nxl_data,nil_data] = size(data{1});
            output(p,q,r) = fhandle(data,nil_data,nxl_data,nt_data);
            %phase(p,q,r) = kpe(data,nil_data,nxl_data,nt_data,t_inc/1000); % Call a function to do something with the data
            fprintf('Completed block %d of %d = %.0f%%\n',r+(nil_block-1)*(q-1)+((nxl_block-1)*(nil_block-1))*(p-1),(nil_block-1)*(nxl_block-1)*(nt_block-1),100*((r+(nil_block-1)*(q-1)+((nxl_block-1)*(nil_block-1))*(p-1))/((nil_block-1)*(nxl_block-1)*(nt_block-1))));
        end
    end
end
end

fclose all;

end

% data = phase;
% ilaxis = round((nil_start+nil_stepout+1:((nil_end-nil_stepout-1)-(nil_stepout+nil_start+1))/(nil_block-2):nil_end-nil_stepout-1));
% xlaxis = round((nxl_start+nxl_stepout+1:((nxl_end-nxl_stepout-1)-(nxl_stepout+nxl_start+1))/(nxl_block-2):nxl_end-nxl_stepout-1));
% taxis = round((nt_start+nt_stepout+1:((nt_end-nt_stepout-1)-(nt_stepout+nt_start+1))/(nt_block-2):nt_end-nt_stepout-1));
% range_in = {ilaxis,xlaxis,taxis};
% range_out = {(1:1:nil),(1:1:nxl),(1:1:nt)};
% 
% clearvars -except data range_in range_out
% 
% [datai] = it2di(phase,range_in,range_out);