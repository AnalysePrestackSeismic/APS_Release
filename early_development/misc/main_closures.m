bitsperbyte = 4;

nil = 295;
nxl = 645;
nt = 838;

contint = 25;

file_in = 'vol_polar_dip.bin';
file_out = 'vol_polar_dip_closure.bin';

vol_in = fread(fopen(file_in),'float32');
vol_in = reshape(vol_in,nt,[])';
vol_out = zeros(nil*nxl,nt);

for nslice = 1:nt
    fprintf('%d of %d slices processing...\n',nslice,nt)
    %slice = slicer(file_in,nslice,nil,nxl,nt,bitsperbyte);
    %slice_closures = findclosures(slice,contint);
    %slice_write(file_out,slice_closures,nslice,nil,nxl,nt,bitsperbyte);
    slice_closures = findclosures(reshape(vol_in(:,nslice),nxl,nil),contint);
    vol_out(:,nslice) = reshape(slice_closures,[],1);
end
    