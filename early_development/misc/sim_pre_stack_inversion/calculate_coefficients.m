function [c1 c2 c3] = calculate_coefficients(angles,no_rho,gamma,k,m)

    c1 = 1+(tan(angles).*tan(angles));
    c2 = (-8)*(gamma^2)*(tan(angles).*tan(angles));
    c3 = ((-0.5)*(tan(angles).*tan(angles)))+(2*(gamma^2)*(sin(angles).*sin(angles)));
    c1 = (0.5.*c1)+(0.5*k.*c2)+(m.*c3);
    c2 = 0.5.*c2;

    if no_rho == 1
        c3=c3*0; % do not invert for density
    end
end