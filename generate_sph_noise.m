function [theta_noise,phi_noise]=generate_sph_noise(theta,phi,mu,sigma)
    psi=normrnd(mu,sigma);
    if psi<0
        psi=-psi;
    end
    % theta_err=psi;
    phi_err=rand*2*pi;
    rot_axis=cross([0,0,1],[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]);
    rot_axis=rot_axis./norm(rot_axis);
    %rot_angle=phi;
    rot_axis_cross=[0,-rot_axis(3),rot_axis(2);-rot_axis(3),0,-rot_axis(1);-rot_axis(2),rot_axis(1),0];
    rot_axis_cross2=rot_axis_cross*rot_axis_cross;
    rot_matrix=eye(3)+sin(theta).*rot_axis_cross+(1-cos(theta)).*rot_axis_cross2;
    err_carte=[sin(psi)*cos(phi_err);sin(psi)*sin(phi_err);cos(psi)];
    rot_err_carte=rot_matrix*err_carte;
    theta_noise=acos(rot_err_carte(3));
    phi_noise=atan(rot_err_carte(2)/rot_err_carte(1));
    if rot_err_carte(1)<0
        phi_noise=phi_noise+pi;
    end
    if phi_noise>pi
        phi_noise=phi_noise-2*pi;
    end
end