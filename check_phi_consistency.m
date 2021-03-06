run_num=1000;
ptsnum=200;
psis=linspace(0,pi/10,100);
% psis=0;
phi_div=36;
% phi_div=1;
chosen_phis=-ones(run_num,length(psis),phi_div,ptsnum+1);
reproj_errs=-ones(run_num,length(psis),phi_div,ptsnum);
chosen_phis_total=cell(run_num,1);
reproj_errs_total=cell(run_num,1);
chosen_Rs=cell(run_num,length(psis),phi_div);
chosen_Ts=cell(run_num,length(psis),phi_div);
epi2_thTooLarge=0;
for m=1:run_num
    [R_gt,T_gt,e1_gt,e2_gt,phi_gt,K,p1OnImage,p2OnImage]=gen_RT;
%     p1_noise=[rand(2,ptsnum).*0.01;zeros(1,ptsnum)];
%     p2_noise=[rand(2,ptsnum).*0.01;zeros(1,ptsnum)];
%     p1_noisy=p1OnImage(:,1:ptsnum)+p1_noise;
%     p2_noisy=p2OnImage(:,1:ptsnum)+p2_noise;
    epi1_sphgt=implane2imsph(e1_gt);
    epi2_sphgt=implane2imsph(e2_gt);
    if epi2_sphgt(1)>pi/2-psis(end)
        epi2_thTooLarge=epi2_thTooLarge+1;
        continue
    end
    rot_matrix=rot2sph_pt(epi2_sphgt);
    for n=1:length(psis)
        psi=psis(n);
        epi2s_sph=gen_sph_iso_pts(rot_matrix,psi,phi_div);
        for o=1:phi_div
            epi2_sph=epi2s_sph(:,o);
            epi12_sph=[epi1_sphgt',epi2_sph'];
            [chosen_R,chosen_T,chosen_phi_avg,chosen_phi]=epipole_corrs_to_RT(epi12_sph,p1OnImage,p2OnImage,ptsnum);
%             [chosen_R,chosen_T,chosen_phi_avg,chosen_phi]=epipole_corrs_to_RT(epi12_sph,p1_noisy,p2_noisy,ptsnum);
            chosen_phis(m,n,o,1)=phi_gt;
            chosen_phis(m,n,o,2:end)=chosen_phi;
            chosen_Rs{m,n,o}=chosen_R;
            chosen_Ts{m,n,o}=chosen_T;
            reproj_errs(m,n,o,:)=compute_reprojection_error_all(chosen_R,chosen_T,p1OnImage,p2OnImage,ptsnum);
%             reproj_errs(m,n,o,:)=compute_reprojection_error_all(chosen_R,chosen_T,p1_noisy,p2_noisy,ptsnum);
        end
    end
    chosen_phis_total{m,1}=chosen_phis;
    reproj_errs_total{m,1}=reproj_errs;
end
% error_runs=[];
% for m=1:run_num
%     maxphi=max(chosen_phis(m,2:end));
%     minphi=min(chosen_phis(m,2:end));
%     if or(maxphi-minphi>1e-9,abs(chosen_phis(m,1)-(maxphi+minphi)/2>1e-9))
%         error_runs(end+1)=m;
%     end
% end
% error_phis=chosen_phis(error_runs,:);