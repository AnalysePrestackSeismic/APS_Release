function [x_dip,y_dip] = calculate_dip(M)
% calculate eigenvalues and eigenvectors
    [V D] = eig(M);
    [l,idx] = sort(abs(diag(D,0)),'descend');
    %l1 = l(1);
    %l2 = l(2);
    %l3 = l(3);
    e1 = V(:,idx(1));
    %e2 = V(:,idx(2));
    %e3 = V(:,idx(3));
    
    x_dip = atand(e1(1)./e1(3));
    y_dip = atand(e1(2)./e1(3));
end