addpath('/data/TZA/dtect/bg_tza_site_survey_vol2_msc/Misc')
file='vol_3D_dip_25_rk_ascii'
data=importdata(file);

siz=size(data)

k=1;

%%%sorting the data takes too long

%for i=1:siz(1)
%    inline=data(i,1);
%    xline=data(i,2);
%    for j=3:siz(2)
%        data_sorted(k,1)=inline;
%       data_sorted(k,2)=xline;
%        data_sorted(k,3)=j-2;
%        data_sorted(k,4)=data(i,j);
%        k=k+1;
%   end
%end
