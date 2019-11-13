% Vary one ground truth epipole, see the spreads of phi and reprojection
% error.
run_num=1e1;
ptsnum=100;
psis=linspace(0,pi/10,101);
epi1sph_max=pi/2;
epi2sph_max=pi/2-pi/10;
% psis=0;
phi_div=10;
% phi_div=1;
chosen_phis=-ones(run_num,length(psis),phi_div,ptsnum+1);
reproj_errs=-ones(run_num,length(psis),phi_div,ptsnum);
chosen_phis_total=cell(run_num,1);
reproj_errs_total=cell(run_num,1);
chosen_Rs=cell(run_num,length(psis),phi_div);
chosen_Ts=cell(run_num,length(psis),phi_div);
epi2_thTooLarge=0;
epi2_thgood=0;
for m=1:run_num
    [R_gt,T_gt,e1_gt,e2_gt,phi_gt,K,p1OnImage,p2OnImage]=gen_RT(epi1sph_max,epi2sph_max);
%     p1_noise=[rand(2,ptsnum).*0.01;zeros(1,ptsnum)];
%     p2_noise=[rand(2,ptsnum).*0.01;zeros(1,ptsnum)];
%     p1_noisy=p1OnImage(:,1:ptsnum)+p1_noise;
%     p2_noisy=p2OnImage(:,1:ptsnum)+p2_noise;
    epi1_sphgt=implane2imsph(e1_gt);
    epi2_sphgt=implane2imsph(e2_gt);
    if epi2_sphgt(1)>epi2sph_max
        epi2_thTooLarge=epi2_thTooLarge+1;
        continue
    else
        epi2_thgood=epi2_thgood+1;
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
            chosen_phis(epi2_thgood,n,o,1)=phi_gt;
            chosen_phis(epi2_thgood,n,o,2:end)=chosen_phi;
            chosen_Rs{epi2_thgood,n,o}=chosen_R;
            chosen_Ts{epi2_thgood,n,o}=chosen_T;
            reproj_errs(epi2_thgood,n,o,:)=compute_reprojection_error_all(chosen_R,chosen_T,p1OnImage,p2OnImage,ptsnum);
%             reproj_errs(m,n,o,:)=compute_reprojection_error_all(chosen_R,chosen_T,p1_noisy,p2_noisy,ptsnum);
        end
    end
end
error_runs=[];
% for m=1:run_num
%     maxphi=max(chosen_phis(m,2:end));
%     minphi=min(chosen_phis(m,2:end));
%     if or(maxphi-minphi>1e-9,abs(chosen_phis(m,1)-(maxphi+minphi)/2>1e-9))
%         error_runs(end+1)=m;
%     end
% end
% error_phis=chosen_phis(error_runs,:);
avg_phistds=-ones(1,length(psis));
avg_phi_diffs=-ones(1,length(psis));
avg_avg_reproj_errs=-ones(1,length(psis));
for m=1:length(psis)
    avg_phistd=0;
    avg_phi_diff=0;
    avg_avg_reproj_err=0;
    chosen_phis_slice=squeeze(chosen_phis(:,m,:,:).*180./pi);
    reproj_errs_slice=squeeze(reproj_errs(:,m,:,:));
    for n=1:run_num
        for o=1:phi_div
            phistd=std(chosen_phis_slice(n,o,2:end));
            avg_phistd=avg_phistd+phistd;
            avg_phi=sum(chosen_phis_slice(n,o,2:end))/ptsnum;
            phi_diff=abs(avg_phi-chosen_phis_slice(n,o,1));
            avg_phi_diff=avg_phi_diff+phi_diff;
            avg_reproj_err=sum(reproj_errs_slice(n,o,:))/ptsnum;
            avg_avg_reproj_err=avg_avg_reproj_err+avg_reproj_err;
        end
    end
    avg_phistd=avg_phistd/run_num/phi_div;
    avg_phistds(m)=avg_phistd;
    avg_phi_diffs(m)=avg_phi_diff/run_num/phi_div;
    avg_avg_reproj_errs(m)=avg_avg_reproj_err.*K(1)/run_num/phi_div;
end
plot(psis.*180./pi,avg_phistds)
xlabel('Epipole perturbation (degree)')
ylabel('Average std(\phi) (degree)')
figure
plot(psis.*180./pi,avg_phi_diffs)
xlabel('Epipole perturbation (degree)')
ylabel('Absolute difference between avg(\phi) and ground truth \phi')
figure
plot(psis.*180./pi,avg_avg_reproj_errs)
xlabel('Epipole perturbation (degree)')
ylabel('Average reprojection error (pixel)')

slices=squeeze(chosen_phis(:,1,1,:));
% slice=slices(:,1,:);
error_runs=[];
for n=1:run_num
    maxphi=max(slices(n,2:end));
    minphi=min(slices(n,2:end));
    if or(maxphi-minphi>1e-9,abs(slices(n,1)-(maxphi+minphi)/2>1e-9))
        error_runs(end+1)=n;
    end
end
error_phis=squeeze(slices(error_runs,:));