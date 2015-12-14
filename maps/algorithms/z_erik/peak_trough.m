inter=[1,2,3,3,2,1,1,2,3,3,2,1,1,2,3]
grad=[4,6,3,7,5,3,9,3,5,7,3,2,6,1,8]
deriv_1_inter=diff(inter)
deriv_2_inter=diff(inter,2)
inter_pt=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
grad_pt=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

for i=1:14
    if deriv_1_inter(i)==0
        if (deriv_2_inter(i)>0|deriv_2_inter(i)<0)
            inter_pt(i)=inter(i);
            grad_pt(i)=grad(i);
        end
    end
end

inter_pt
grad_pt
