inter=[1,2,3,3,2,1,1,2,3,3,2,1,1,2,3]
grad=[4,6,3,7,5,3,9,3,5,7,3,2,6,1,8]
deriv_1_inter=diff(inter)
inter_t=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
grad_t=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
sign_d1=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

for i=2:14
    if deriv_1_inter(i)>0
        sign_d1(i)=1;
    elseif deriv_1_inter(i)<0
        sign_d1(i)=-1;
    end
    if sign_d1(i-1)<sign_d1(i)
            inter_t(i)=inter(i);
			grad_t(i)=grad(i);
    end
end

sign_d1
inter_t
grad_t