% traces2cube

n_iline=segy{1}.n_iline;
n_xline=segy{1}.n_xline;
n_samples=500;  %segy{1}.n_samples;
data=traces.data;

min_iline=1;
max_iline=n_iline;
min_xline=1;
max_xline=n_xline;

cube=zeros(n_samples,max_iline-min_iline+1,max_xline-min_xline+1);

for i=min_iline:max_iline
    i
    ind=(min_xline+(i-1)*n_xline):1:(max_xline+(i-1)*n_xline);
    cube(:,i-min_iline+1,:)=data(1:n_samples,ind);
end