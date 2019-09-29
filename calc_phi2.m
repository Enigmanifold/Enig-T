function [phi1,phi2,phi1_alt,phi2_alt]=calc_phi2(p1,p2,e1,e2,ptsnum)
    gamma_cross_e=zeros(3,ptsnum);
    e_norm=norm(e1);
    ebar_norm=norm(e2);
    ebar_dot_gamma=zeros(1,ptsnum);
    e_dot_ebar=e1'*e2;
    eebar_norm=e_norm*ebar_norm;
    eebar_ebare=e_norm.*e2+ebar_norm.*e1;
    ebar_cross_e=cross(e2,e1);
    eeee=eebar_norm+e_dot_ebar;
    eeeee=ebar_norm*eeee;
    R2_additional_axis_y=-(e2(1)+1)/e2(2);
    R2_additional_cross=[0,-1,R2_additional_axis_y;1,0,-1;-R2_additional_axis_y,1,0];
    R2_additional_cross2=R2_additional_cross*R2_additional_cross;
    R2_additional=eye(3)+2*e2(2)^2/(2*e2(2)^2+e2(1)^2+2*e2(1)+1)*R2_additional_cross2;
    phi1=zeros(1,ptsnum);
    phi1_alt=zeros(1,ptsnum);
    for m=1:ptsnum
        gamma_cross_e(:,m)=cross(p1(:,m),e1);
        ebar_dot_gamma(m)=e2'*p1(:,m);
    end
    for m=1:ptsnum
        tanphi1=-(p2(:,m)'*(eeeee.*gamma_cross_e(:,m)+p1(:,m)'*ebar_cross_e.*eebar_ebare))/(p2(:,m)'*(eeee.*(e_dot_ebar.*p1(:,m)-ebar_dot_gamma(m).*e1)+p1(:,m)'*ebar_cross_e.*ebar_cross_e));
        ph1_new=atan(tanphi1);
        if ph1_new<0
            ph1_new=ph1_new+pi;
        end
        phi1(m)=ph1_new;
    end
end