function [data,cube,Lcube_true]=model4(xmax,ymax,zmax,nfaults,amax,bmax,noise_frac)

[data,cube,Lcube_true]=model3(xmax,ymax,zmax,nfaults,amax,bmax);

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