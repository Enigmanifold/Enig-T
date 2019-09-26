run('epipole_theta_to_RT.m');
gamma_cross_e=zeros(3,ptsnum);
gamma_cross_ebar=zeros(3,ptsnum);
e_coeff=zeros(1,ptsnum);
e_norm=norm(e1);
ebar_norm=norm(e2);
e_dot_ebar=e1'*e2;
phi1_new=zeros(1,ptsnum);
phi1_new_alt=zeros(1,ptsnum);
eebar_norm=e_norm*ebar_norm;
R2_additional_axis_y=-(e2(1)+1)/e2(2);
R2_additional_cross=[0,-1,R2_additional_axis_y;1,0,-1;-R2_additional_axis_y,1,0];
R2_additional_cross2=R2_additional_cross*R2_additional_cross;
R2_additional=eye(3)+2*e2(2)^2/(2*e2(2)^2+e2(1)^2+2*e2(1)+1)*R2_additional_cross2;
for m=1:ptsnum
    gamma_cross_e(:,m)=cross(p1(:,m),e1);
    gamma_cross_ebar(:,m)=cross(p1(:,m),e2);
    e_coeff(m)=e1'*gamma_cross_ebar(:,m)./eebar_norm;
end
e_norm2=e_norm^2;
e_perp_proj=e2-e_dot_ebar./e_norm2.*e1;
epi_fiesta=e_norm2/(eebar_norm+e_dot_ebar);
tanphi1=zeros(1,ptsnum);
for m=1:ptsnum
    N2_unnormalized=gamma_cross_e(:,m)+e_coeff(m).*e1+epi_fiesta.*e_coeff(m).*e_perp_proj;
    numer=-p2_to_e2(:,m)'*N2_unnormalized*ebar_norm;
    denom=p2_to_e2(:,m)'*cross(e2,N2_unnormalized);
    tanphi1=numer/denom;
    ph1_new=atan(tanphi1);
    if ph1_new<0
        ph1_new=ph1_new+pi;
    end
    phi1_new(m)=ph1_new;    
    N2_unnormalized_alt=R2_additional*N2_unnormalized;
    numer_alt=-p2_to_e2(:,m)'*N2_unnormalized_alt*ebar_norm;
    denom_alt=p2_to_e2(:,m)'*cross(e2,N2_unnormalized_alt);
    tanphi1_alt=numer_alt/denom_alt;
    ph1_new_alt=atan(tanphi1_alt);
    if ph1_new_alt<0
        ph1_new_alt=ph1_new_alt+pi;
    end
    phi1_new_alt(m)=ph1_new_alt;
end
phi2_new=phi1_new+repmat(pi,1,ptsnum);
phi2_new_alt=phi1_new_alt+repmat(pi,1,ptsnum);