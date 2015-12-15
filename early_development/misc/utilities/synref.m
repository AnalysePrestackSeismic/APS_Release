function [ref, imp] = synref(nt,sparsity,maxamp)
    ref = zeros(nt,1);
    imp = zeros(nt,1);
    imp(1) = 1500;
    while or(min(imp)<1500,max(imp)>15000)
        for ii=1:nt
            if rand > sparsity
                ref(ii) = randn;
            else
                ref(ii) = 0;
            end
        end
        ref = maxamp*(ref/(max(abs(ref))));
        for jj=1:nt-1
            imp(jj+1) = imp(jj)*((1+ref(jj))/(1-ref(jj)));
        end
    end
end