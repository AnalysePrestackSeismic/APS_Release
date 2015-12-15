function param_store=fit_2d_morlet_test2(t,w,vals,search,res_min,num_mors,inline,xline)


min_res=res_min+1;


n=1;
new_res=vals/255;
total_res=sum(sum(new_res));

while n<=num_mors

    maxi=max(max(new_res));
    [r,c]=find(new_res==(max(max(new_res))));
    if n==1; fir_maxi=maxi;end
    
    old_res=min_res;
    
    if maxi/fir_maxi<0.04;break;end%disp('No maximum greater than 5% of the original maximum found, exiting loops...');break;end
    %disp(min_res)
   
    marker=0;
    for i=(r-1):-1:2
        if (abs(vals(i,c)-maxi/2)<abs(vals(i-1,c)-maxi/2)) & (abs(vals(i,c)-maxi/2)<abs(vals(i+1,c)-maxi/2))
            tsig(1)=abs(i-r);
            marker=marker+1;break;
        end
        if marker==0;tsig(1)=abs(1-r);end
    end
    marker=0;
   for i=(r+1):(size(vals,1)-1)
        if (abs(vals(i,c)-maxi/exp(2))<abs(vals(i-1,c)-maxi/2)) & (abs(vals(i,c)-maxi/2)<abs(vals(i+1,c)-maxi/2))
            tsig(2)=abs(i-r);
            marker=marker+1;break;
        end
        if marker==0;tsig(2)=abs(r-size(vals,1));end
    end    
    marker=0;
    for i=(c-1):-1:2
        if (abs(vals(r,i)-maxi/2)<abs(vals(r,i-1)-maxi/2)) & (abs(vals(r,i)-maxi/2)<abs(vals(r,i+1)-maxi/2))
            fsig(1)=abs(i-c);
            marker=marker+1;break;
        end
        if marker==0;fsig(1)=abs(c-1);end
    end
    marker=0;
    for i=(c+1):(size(vals,2)-1)
        if (abs(vals(r,i)-maxi/2)<abs(vals(r,i-1)-maxi/2)) & (abs(vals(r,i)-maxi/2)<abs(vals(r,i+1)-maxi/2))
            fsig(2)=abs(i-c);
            marker=marker+1;break;
        end
        if marker==0;fsig(2)=abs(c-size(vals,2));end
    end
    
    fsig=fsig*(w(2)-w(1));
    tsig=tsig*(t(2)-t(1));

        sigt=sqrt((-1*min(tsig)^2)/(14*log(0.5)));   %%% for synthetic
     sigw=sqrt((-1*min(fsig)^2)/(16*log(0.5)));   %%% for synthetic



% sigt=9.5;
% sigw=8.4;

    lowest_res_params=[t(r),w(c),sigt,sigw,maxi];
   
    
while min_res>res_min
m=1;

   u=t(r);
   wm=w(c);
   amp=maxi;
           for ksigt=1:3 %%%% sigt
            sigt=lowest_res_params(3)+(ksigt-2)*search.sigt;
                for ksigw=1:3 %%%% sigw
                    sigw=lowest_res_params(4)+(ksigw-2)*search.sigw;
                    
                    res(m)=0;
                     for i=1:size(t,2)
                        for j=1:size(w,2)

                            gauss2=amp*exp(-((((w(j)-wm)^2)/(2*sigw^2))+(((t(i)-u)^2)/(2*sigt^2))));
                            if gauss2<1e-28;gauss2=0;end
                            res(m)=res(m)+(abs(new_res(i,j)-gauss2));
                        end
                     end
                 
                    param(m,1:5)=[u,wm,sigt,sigw,amp];
                    m=m+1;
                end
           end

[min_res,loc_min_res]=min(res);
min_res;

if lowest_res_params==param(loc_min_res,1:5);break;end%disp('Found minimum, exiting loop...');break;end
if (abs((old_res-min_res))/old_res)<0.005;break;end%disp('Residual not decreasing, exiting loop...');break;end

lowest_res_params=param(loc_min_res,1:5);

end

for i=1:size(t,2)
    for j=1:size(w,2)

        gauss2=lowest_res_params(5)*exp(-((((w(j)-lowest_res_params(2)).^2)./(2*lowest_res_params(4)^2))+(((t(i)-lowest_res_params(1)).^2)./(2*lowest_res_params(3)^2))));
        new_res(i,j)=new_res(i,j)-gauss2;

    end
end                    


% figure(n);
% a=axes;
% set(a,'FontSize',18);
% %imagesc(w,t,new_res,[-0.2 0.05])
% imagesc(w,t,new_res)
% xlabel('Frequency (Hz)','FontSize',18)
% ylabel('Time (msec)','FontSize',18)
% c=colorbar;
% set(c,'FontSize',18);
% ylabel(c,'Amplitude','FontSize',18)


param_store(n,:)=[lowest_res_params(1), lowest_res_params(2),lowest_res_params(3),lowest_res_params(4),lowest_res_params(5),inline,xline,min_res];

n=n+1;

end

end


