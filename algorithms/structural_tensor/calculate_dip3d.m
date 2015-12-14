function [nx,ny,nz,planarity] = calculate_dip3d(Ix,Iy,Iz,sigma,scale_sigma)
% calculate eigenvalues and eigenvector
    Ixy = Ix.*Iy;
    Ixz = Ix.*Iz;
    Iyz = Iy.*Iz;
    Iyy = Iy.*Iy;
    Ixx = Ix.*Ix;
    Izz = Iz.*Iz;

    % Smooth gradient components with a Gaussian
    Ixy = imgaussian(Ixy,sigma,scale_sigma*sigma);
    Ixz = imgaussian(Ixz,sigma,scale_sigma*sigma);
    Iyz = imgaussian(Iyz,sigma,scale_sigma*sigma);
    Iyy = imgaussian(Iyy,sigma,scale_sigma*sigma);
    Ixx = imgaussian(Ixx,sigma,scale_sigma*sigma);
    Izz = imgaussian(Izz,sigma,scale_sigma*sigma);

    Ixy = Ixy(:);
    Ixz = Ixz(:);
    Iyz = Iyz(:);
    Iyy = Iyy(:);
    Ixx = Ixx(:);
    Izz = Izz(:);
    
    parfor i_loop = 1:numel(Ixx);
        [V,D] = eig([Ixx(i_loop), Ixy(i_loop) , Ixz(i_loop); Ixy(i_loop) , Iyy(i_loop) , Iyz(i_loop); Ixz(i_loop) , Iyz(i_loop) , Izz(i_loop)]);   
        [l,idx] = sort(abs(diag(D,0)),'descend');
        l1 = l(1);
        l2 = l(2);
        %l3 = l(3);
        e1 = V(:,idx(1));
        %e2 = V(:,idx(2));
        %e3 = V(:,idx(3));
        
        planarity(i_loop,1) = (l1-l2)./l1; % close to 1 where events are locally planar
        nx(i_loop,1) = e1(1);
        ny(i_loop,1) = e1(2);
        nz(i_loop,1) = e1(3);        
    end
end