function [ref] = reflectivity(t_ref, x_ref, sparsity)

ref=zeros(t_ref:x_ref);
for i=1:t_ref
    if sparsity<rand
        ref(i,:)=randn*0.1;
    else
        ref(i,:)=0;
    end
end

