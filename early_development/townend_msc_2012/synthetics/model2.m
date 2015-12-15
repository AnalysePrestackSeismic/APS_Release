function [data,cube,Lcube_true]=model2(xmax,ymax,zmax,xc,yc,zc,a,b,theta,phi,PHI,noise_frac)

[data,cube,Lcube_true]=model1(xmax,ymax,zmax,xc,yc,zc,a,b,theta,phi,PHI);

for trace=1:size(data,2);
    for z=1:size(data,1);
        rnd=rand(1);
        if rnd<=noise_frac && data(z,trace)==10
            data(z,trace)=0;
        elseif rnd<=noise_frac && data(z,trace)~=10
            data(z,trace)=10;
        end
    end
end