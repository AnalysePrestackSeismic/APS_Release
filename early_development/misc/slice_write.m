function [] = slice_write(file_out,slice_closures,nslice,nil,nxl,nt,bitsperbyte)

global fid2;
slice_closures = reshape(slice_closures,nil*nxl,[]);
slice_closures = full(slice_closures);

if nslice == 1
    fid2 = fopen(file_out,'w');
end
fseek(fid2,(nslice-1)*bitsperbyte,'bof');
for i = 1:nil*nxl
    fwrite(fid2,slice_closures(i,1),'float32',bitsperbyte*(nt-1));
end
if nslice == nt
    fclose(fid2);
end
end

