function [slice] = slicer(file_in,nslice,nil,nxl,nt,bitsperbyte)

slice = zeros(nil*nxl,1);
fid1 = fopen(file_in);
fseek(fid1,(nslice-1)*bitsperbyte,'bof');
for i = 1:nil*nxl
    slice(i,1) = fread(fid1,1,'float32',bitsperbyte*(nt-1));
end
fclose(fid1);
slice = reshape(slice,nxl,[]);
end

