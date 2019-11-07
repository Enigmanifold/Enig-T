run_num=1000;
ptsnum=100;
% psis=linspace(0,pi/20,100);
% psis=0;
% phi_div=36;
noise_levels=linspace(0,10,101);
% phi_div=1;
chosen_phis=-ones(length(noise_levels),run_num,ptsnum+1);
reproj_errs=-ones(length(noise_levels),run_num,ptsnum);
% chosen_phis_total=cell(length(noise_levels),1);
% reproj_errs_total=cell(length(noise_levels),1);
% chosen_Rs=cell(run_num,length(psis),phi_div);
% chosen_Ts=cell(run_num,length(psis),phi_div);
% epi2_thTooLarge=0;
for m=1:length(noise_levels)
    noise_level=noise_levels(m);
    for n=1:run_num
        [R_gt,T_gt,e1_gt,e2_gt,phi_gt,K,p1,p2,p1_pixel,p2_pixel]=gen_RT;
%         invK=inv(K);
        randr=normrnd(0,noise_level,[2,ptsnum]); % Introduce noise.
        randtheta=rand(2,ptsnum).*2.*pi;
        p1_noise=[randr(1,:).*cos(randtheta(1,:));randr(1,:).*sin(randtheta(1,:));zeros(1,ptsnum)];
        p2_noise=[randr(2,:).*cos(randtheta(2,:));randr(2,:).*sin(randtheta(2,:));zeros(1,ptsnum)];
        p1_noisy_pix=p1_pixel(:,1:ptsnum)+p1_noise;
        p2_noisy_pix=p2_pixel(:,1:ptsnum)+p2_noise;
        p1_noisy_fl=K\p1_noisy_pix;
        p2_noisy_fl=K\p2_noisy_pix;
        epi1_sphgt=implane2imsph(e1_gt);
        epi2_sphgt=implane2imsph(e2_gt);
        epi12_sph=[epi1_sphgt',epi2_sphgt'];
        [chosen_R,chosen_T,chosen_phi_avg,chosen_phi]=epipole_corrs_to_RT(epi12_sph,p1_noisy_fl,p2_noisy_fl,ptsnum);
    %     [chosen_R,chosen_T,chosen_phi_avg,chosen_phi]=epipole_corrs_to_RT(epi12_sph,p1_noisy,p2_noisy,ptsnum);
        chosen_phis(m,n,1)=phi_gt;
        chosen_phis(m,n,2:end)=chosen_phi;
%         chosen_Rs{n,1}=chosen_R;
%         chosen_Ts{n,1}=chosen_T;
        reproj_errs(m,n,:)=compute_reprojection_error_all(chosen_R,chosen_T,p1_noisy_pix,p2_noisy_pix,ptsnum,K);
    %     reproj_errs(m,n,o,:)=compute_reprojection_error_all(chosen_R,chosen_T,p1_noisy,p2_noisy,ptsnum);
    end
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