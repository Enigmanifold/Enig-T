run_num=100000;
ptsnum=20;
psis=linspace(0,pi/20,100);
% psis=0;
phi_div=36;
% phi_div=1;
chosen_phis=-ones(run_num,ptsnum+1);
reproj_errs=-ones(run_num,ptsnum);
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
    epi12_sph=[epi1_sphgt',epi2_sphgt'];
    [chosen_R,chosen_T,chosen_phi_avg,chosen_phi]=epipole_corrs_to_RT(epi12_sph,p1OnImage,p2OnImage,ptsnum);
%     [chosen_R,chosen_T,chosen_phi_avg,chosen_phi]=epipole_corrs_to_RT(epi12_sph,p1_noisy,p2_noisy,ptsnum);
    chosen_phis(m,1)=phi_gt;
    chosen_phis(m,2:end)=chosen_phi;
    chosen_Rs{m,1}=chosen_R;
    chosen_Ts{m,1}=chosen_T;
    reproj_errs(m,:)=compute_reprojection_error_all(chosen_R,chosen_T,p1OnImage,p2OnImage,ptsnum);
%     reproj_errs(m,n,o,:)=compute_reprojection_error_all(chosen_R,chosen_T,p1_noisy,p2_noisy,ptsnum);
end
error_runs=[];
for m=1:run_num
    maxphi=max(chosen_phis(m,2:end));
    minphi=min(chosen_phis(m,2:end));
    if or(maxphi-minphi>1e-9,abs(chosen_phis(m,1)-(maxphi+minphi)/2>1e-9))
        error_runs(end+1)=m;
    end
end
error_phis=chosen_phis(error_runs,:);