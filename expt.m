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









run('setup_episearch.m');
% results=episearch(epi12_sph,p1,p2,ptsnum,top_frac,method,measure,diff_amplifier,options);
top_frac=0.1;
ptsnum=30;
diff_amplifier=100000;
options_GDstdphi=optimset('MaxFunEvals',10000,'MaxIter',100000,'TolFun',1e-4);
func=@(x)calc_phi3_multi(x,p1,p2,ptsnum,diff_amplifier);