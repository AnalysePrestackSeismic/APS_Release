%function [filename,filename_out,md_rs,age_rs] = LogInterpolation(well,md,age)
%LOGINTERPOLATION Summary of this function goes here
%   Detailed explanation goes here

filename='Chaza1_md_age.txt';
filename_out='Chaza1_md_age_interpolated.txt';
[well md age]=textread(filename,'%s %f %f');
type='linear';

min_val=min(md);
max_val=max(md);

delta=0.1524;

md_rs=interp1(md,md,min_val:delta:max_val,type)';
age_rs=interp1(md,age,min_val:delta:max_val,type)';
out=[md_rs age_rs];
dlmwrite(filename_out,out,'delimiter',' ','precision',8);

%end

