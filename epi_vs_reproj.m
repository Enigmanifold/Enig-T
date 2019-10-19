% Generate random pose, using GT phi and candidate epipoles to calculate
% reprojection errors.
run('epipole_theta_to_RT.m');
if ~isequal(min(ind),max(ind))
    disp('inconsistent phi calculation.')
    return
end
switch min(ind)
    case 1
        phi=phi1(1);
    case 2
        phi=phi2(1);
    case 3
        phi=phi1_alt(1);
    case 4
        phi=phi2_alt(1);
end
[epi1_sph(1),epi1_sph(2)]=implane2imsphere(e1);
[epi2_sph(1),epi2_sph(2)]=implane2imsphere(e2);
maxpsi=pi/2-epi2_sph(1);
if maxpsi<0.5*pi/2
    disp('maxpsi too small.');
    return
end
alpha_div=100;
beta_div=501;
alphas1=linspace(0,2*pi*(1-1/alpha_div),alpha_div);
betas1=linspace(0,0.5*pi/2,beta_div);

epi2_sph_candis=cell(beta_div);
reproj_errs={};
reproj_errs_array=-ones(beta_div,alpha_div);
for m=1:beta_div
   psi1=betas1(m);
   epi2_sph_candi=gen_sph_iso_pts(epi2_sph,psi1,alpha_div);
   epi2_sph_candis{m}=epi2_sph_candi;
end
for m=1:beta_div
    epi2_sph_ring=epi2_sph_candis{m};
    for n=1:alpha_div
        epi2_sph_pt=epi2_sph_ring(:,n);
        reproj_err=compute_reprojection_error(epi1_sph,epi2_sph_pt,phi,p1(:,1:30),p2(:,1:30));
        reproj_errs{end+1,1}=epi2_sph_pt;
        reproj_errs{end,2}=betas1(m);
        reproj_errs{end,3}=reproj_err;
        reproj_errs{end,4}=sum(reproj_err)./length(reproj_err);
        reproj_errs_array(m,n)=reproj_errs{end,4};
    end
end
% reproj_errs_array(:,1)=repmat(linspace(0,2*pi*(1-1/alpha_div),alpha_div)',beta_div,1);
% for m=1:beta_div*alpha_div
%     reproj_errs_array(m,2)=reproj_errs{m,2};
%     reproj_errs_array(m,3)=reproj_errs{m,4};
% end

[Alpha,Beta]=meshgrid(alphas1,betas1);
sf=surf(Alpha.*180./pi,Beta.*180./pi,reproj_errs_array);
% sf.EdgeColor='none';
title('Reprojection error w.r.t. (\alpha,\beta)')
xlabel('\alpha (degree)')
ylabel('\beta (degree)')
zlabel('s (focal length)')
% ang_dists=-ones(psi_div,phi_div);
% for m=1:psi_div
%     epi2_sph_ring=epi2_sph_candis{m};
%     for n=1:phi_div
%         ang_dists(m,n)=calculate_ang_dist(epi2_sph,epi2_sph_ring(:,n));
%     end
% end
% figure
% for m=1:phi_div
%     epi2_pt=epi2_sph_ring(:,m);
%     scatter3(sin(epi2_pt(1))*cos(epi2_pt(2)),sin(epi2_pt(1))*sin(epi2_pt(2)),cos(epi2_pt(1)));
%     hold on
% end
% hold on
% scatter3(sin(epi2_sph(1))*cos(epi2_sph(2)),sin(epi2_sph(1))*sin(epi2_sph(2)),cos(epi2_sph(1)));
% sph_pt=epi2_sph;
% rot_axis=cross([0,0,1],[sin(sph_pt(1))*cos(sph_pt(2)),sin(sph_pt(1))*sin(sph_pt(2)),cos(sph_pt(1))]);
% rot_axis=rot_axis./norm(rot_axis);
% %rot_angle=phi;
% rot_axis_cross=cross_matrix(rot_axis);
% rot_axis_cross2=rot_axis_cross*rot_axis_cross;
% rot_matrix=eye(3)+sin(sph_pt(1)).*rot_axis_cross+(1-cos(sph_pt(1))).*rot_axis_cross2;
% rot001=rot_matrix*[0;0;1];
% rot001=rot001./rot001(3);
% inv_epi2_carte_ring=zeros(3,phi_div);
% inv_epi2_sph_ring=zeros(2,phi_div);
% inv_ang_dist=-ones(1,phi_div);
% for m=1:phi_div
%     inv_epi2_carte_ring(:,m)=inv(rot_matrix)*[sin(epi2_sph_ring(1,m))*cos(epi2_sph_ring(2,m)),sin(epi2_sph_ring(1,m))*sin(epi2_sph_ring(2,m)),cos(epi2_sph_ring(1,m))]';
%     vec=inv_epi2_carte_ring(:,m)./inv_epi2_carte_ring(3,m);
%     [inv_epi2_sph(1,m),inv_epi2_sph(2,m)]=implane2imsphere(vec);
%     inv_ang_dist(m)=calculate_ang_dist([0;0],inv_epi2_sph(:,m));
% end
% figure
% scatter3(inv_epi2_carte_ring(1,:),inv_epi2_carte_ring(2,:),inv_epi2_carte_ring(3,:))
% hold on
% scatter3(0,0,1)