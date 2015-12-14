function [g] = tensfiltA(data,nX,nZ,u2)
% apply A'A filter from D. Hale, 2007, CWP Report 567, "Local dip filtering
% with directional Laplacians"
    for ii = 2:nZ
        for jj=2:nX
            u2i = u2(ii,jj);
            u1i = sqrt(1-u2i*u2i);
            a11 = 1-u1i*u1i;
            a12 = -u1i*u2i;
            a22 = 1-u2i*u2i;
            fa = data(ii,jj)-data(ii-1,jj-1);
            fb = data(ii,jj-1)-data(ii-1,jj);
            f1 = 0.5*(fa-fb);
            f2 = 0.5*(fa+fb);
            g1 = a11*f1+a12*f2;
            g2 = a12*f1+a22*f2;
            ga = 0.5*(g1+g2);
            gb = 0.5*(g1-g2);
            g(ii,jj) = ga;
            g(ii-1,jj-1) = g(ii-1,jj-1)-ga;
            g(ii,jj-1) = g(ii,jj-1)-gb;
            g(ii-1,jj) = g(ii-1,jj)+gb;
        end
    end
end