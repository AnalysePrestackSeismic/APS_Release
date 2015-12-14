function [p] = calculate_dip(Ix,Iz,sigma,scale_sigma)
% calculate eigenvalues and eigenvector
    Ixz = Ix.*Iz;
    Ixx = Ix.*Ix;
    Izz = Iz.*Iz;
    
    Ixz = imgaussian(Ixz,sigma,scale_sigma*sigma);
    Ixx = imgaussian(Ixx,sigma,scale_sigma*sigma);
    Izz = imgaussian(Izz,sigma,scale_sigma*sigma);
    

    parfor i_loop = 1:numel(Ixx);
        [V D] = eig([Ixx(i_loop), Ixz(i_loop); Ixz(i_loop), Izz(i_loop)]);
        [l,idx] = sort(abs(diag(D,0)),'descend');
        %l1 = l(1);
        %l2 = l(2);
        %l3 = l(3);
        e1 = V(:,idx(1));
        %e2 = V(:,idx(2));
        %e3 = V(:,idx(3));

        p(i_loop) = atand(e1(1)/e1(2));
        %q(i_loop) = atand(e1(2)/e1(3));
    end
end