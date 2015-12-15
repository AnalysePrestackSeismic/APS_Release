fhandle = @kpe;

filenames{1,1} = '/seismic/MAD/dtect_projects/rsimm_test/near.bin';

bitsperbyte = 4;

nil = 101;
nxl = 101;
nt = 626;

nil_stepout = 50;
nxl_stepout = 50;
nt_stepout = 312;

il_olap = 1;
xl_olap = 1;
t_olap = 1;

[pvol] = minivolread(fhandle,filenames,bitsperbyte,nil,nxl,nt,nil_stepout,nxl_stepout,nt_stepout,il_olap,xl_olap,t_olap);