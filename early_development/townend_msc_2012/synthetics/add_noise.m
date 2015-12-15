%function [data,cube,Lcube_true]=model2(xmax,ymax,zmax,xc,yc,zc,a,b,theta,phi,PHI,noise_frac)

%[data,cube,Lcube_true]=model1(xmax,ymax,zmax,xc,yc,zc,a,b,theta,phi,PHI);

noise_frac=0.2;

for z=1:size(Lcube_true,1);
    for x=1:size(Lcube_true,2);
        for y=1:size(Lcube_true,3);
            rnd=rand(1);
            if rnd<=noise_frac && Lcube_true(z,x,y)>0
                Lcube_true(z,x,y)=0;
            elseif rnd<=noise_frac && Lcube_true(z,x,y)==0
                Lcube_true(z,x,y)=10;
            end
        end
    end
end