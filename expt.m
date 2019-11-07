% epi2_sph_pts=zeros(2,beta_div*alpha_div);
% epi2_carte_pts=zeros(3,beta_div*alpha_div);
% for m=1:beta_div*alpha_div
%     epi2_sph_pts(:,m)=reproj_errs{m,1};
%     epi2_carte_pts(:,m)=[sin(epi2_sph_pts(1,m))*cos(epi2_sph_pts(2,m));sin(epi2_sph_pts(1,m))*sin(epi2_sph_pts(2,m));cos(epi2_sph_pts(1,m))];
% end
% for m=1:beta_div
%     scatter3(epi2_carte_pts(1,(m-1)*alpha_div+1:m*alpha_div),epi2_carte_pts(2,(m-1)*alpha_div+1:m*alpha_div),epi2_carte_pts(3,(m-1)*alpha_div+1:m*alpha_div),'filled')
%     hold on
% end
% scatter3(v2(1),v2(2),v2(3),'*','k')
% title('isometric epipole sampling on image sphere')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% % legend('angular distance 0(coincides with ground truth)','angular distance 0.05pi','angular distance 0.1pi','angular distance 0.15pi','angular distance 0.2pi','angular distance 0.25pi','ground truth')
% hold off


% run('generate_RT.m');
% [epi1_sph(1),epi1_sph(2)]=implane2imsphere(e1);
% [epi2_sph(1),epi2_sph(2)]=implane2imsphere(e2);
% epi12_sph=[epi1_sph,epi2_sph];
% ptsnum=30;
% [chosen_R,chosen_T,chosen_phi]=epipole_corrs_to_RT(epi12_sph,p1,p2,ptsnum);
% invR=inv(R);
% R_dists=-ones(1,1000);
% for m=1:10000
%     R_dists(m)=acos((trace(invR*chosen_R)-1)/2);
% end
% max_R_dist=max(R_dists);


lb=0;
ub=2*pi;
lb2=[0,-pi,0,-pi,0];
ub2=[pi/2-0.0001,pi,pi/2-0.001,pi,2*pi];
nvars=1;
nvars2=5;
options=optimoptions('particleswarm','FunctionTolerance',1e-12,'MaxStallIterations',100,'SwarmSize',100);
options2=optimoptions('particleswarm','FunctionTolerance',1e-12,'MaxStallIterations',200,'SwarmSize',3000);
func=@(x)norm(phi_to_RE(x,epi12_sphgt,p1_pixel,p2_pixel,ptsnum,K));
[best_phi,calc_RE_gt]=particleswarm(func,nvars,lb,ub,options);
func2=@(x)norm(epiphi_to_RE(x,p1_pixel,p2_pixel,ptsnum,K));
% [best_epiphi,calc_RE_gt2]=particleswarm(func2,nvars2,lb2,ub2,options2);
e1_searched=-ones(2,1);
e2_searched=-ones(2,1);
[e1_searched(1),e1_searched(2)]=imsphere2implane(results_epiphi4(1),results_epiphi4(2));
[e2_searched(1),e2_searched(2)]=imsphere2implane(results_epiphi4(3),results_epiphi4(4));
phi_searched=results_epiphi4(5);
% [e1_searched(1),e1_searched(2)]=imsphere2implane(epiphi_gt(1),epiphi_gt(2));
% [e2_searched(1),e2_searched(2)]=imsphere2implane(epiphi_gt(3),epiphi_gt(4));
% phi_searched=epiphi_gt(5);
e1_searched=[e1_searched;1];
e2_searched=[e2_searched;1];

% [R_searched,R_alt_searched,T_searched,T_alt_searched]=epipoles_phi_to_RT(e1,e2,best_phi,1);
[R_searched,R_alt_searched,T_searched,T_alt_searched]=epipoles_phi_to_RT(e1_searched,e2_searched,phi_searched,1);
RE_searched1=compute_reprojection_error_all(R_searched,T_searched,p1_pixel_all(:,valid_pairs),p2_pixel_all(:,valid_pairs),ptsnum,K);
RE_searched2=compute_reprojection_error_all(R_alt_searched,T_searched,p1_pixel_all(:,valid_pairs),p2_pixel_all(:,valid_pairs),ptsnum,K);
RE_searched3=compute_reprojection_error_all(R_searched,T_alt_searched,p1_pixel_all(:,valid_pairs),p2_pixel_all(:,valid_pairs),ptsnum,K);
RE_searched4=compute_reprojection_error_all(R_alt_searched,T_alt_searched,p1_pixel_all(:,valid_pairs),p2_pixel_all(:,valid_pairs),ptsnum,K);
RE_searched1_norm=norm(RE_searched1);
RE_searched2_norm=norm(RE_searched2);
RE_searched3_norm=norm(RE_searched3);
RE_searched4_norm=norm(RE_searched4);
[RE_final,RE_ind]=min([RE_searched1_norm,RE_searched2_norm,RE_searched3_norm,RE_searched4_norm]);
switch RE_ind
    case 1
        R_final=R_searched;
        T_final=T_searched;
    case 2
        R_final=R_alt_searched;
        T_final=T_searched;
    case 3
        R_final=R_searched;
        T_final=T_alt_searched;
    case 4
        R_final=R_alt_searched;
        T_final=T_alt_searched;
end
% [best_R,best_T]=epipoles_phi_to_RT(epi12_sphgt,