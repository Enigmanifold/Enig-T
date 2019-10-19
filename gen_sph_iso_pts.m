function sph_pts=gen_sph_iso_pts(rot_matrix,psi,ptsnum)
% Generate spherically isometric points w.r.t. theta=0, then rotate
% according to rot_matrix.
    sph_pts=zeros(2,ptsnum);
    % theta_err=psi;
    phis=linspace(0,2*pi*(1-1/ptsnum),ptsnum);
%     rot_axis=cross([0,0,1],[sin(sph_pt(1))*cos(sph_pt(2)),sin(sph_pt(1))*sin(sph_pt(2)),cos(sph_pt(1))]);
%     rot_axis=rot_axis./norm(rot_axis);
%     %rot_angle=phi;
%     rot_axis_cross=cross_matrix(rot_axis);
%     rot_axis_cross2=rot_axis_cross*rot_axis_cross;
%     rot_matrix=eye(3)+sin(sph_pt(1)).*rot_axis_cross+(1-cos(sph_pt(1))).*rot_axis_cross2;
    for m=1:ptsnum
        dev_carte=[sin(psi)*cos(phis(m));sin(psi)*sin(phis(m));cos(psi)];
        rot_dev_carte=rot_matrix*dev_carte;
        rot_dev_carte=rot_dev_carte./norm(rot_dev_carte);
        theta_dev=acos(rot_dev_carte(3));
        phi_dev=atan(rot_dev_carte(2)/rot_dev_carte(1));
        if rot_dev_carte(1)<0
            phi_dev=phi_dev+pi;
        end
        if phi_dev>pi
            phi_dev=phi_dev-2*pi;
        end
        sph_pts(:,m)=[theta_dev;phi_dev];
    end
end