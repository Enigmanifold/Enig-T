function min_avg_reproj_err=epipole_corrs_to_RE(epi12_sph,p1,p2,ptsnum,sm)
    [~,~,chosen_phi]=epipole_corrs_to_RT(epi12_sph,p1,p2,ptsnum);
    epi1_sph=epi12_sph(1:2);
    epi2_sph=epi12_sph(3:4);
    p1_samples=p1(:,1:ptsnum);
    p2_samples=p2(:,1:ptsnum);
    if sm=='s'
        min_avg_reproj_err=compute_reprojection_error_avg(chosen_phi,epi1_sph,epi2_sph,p1_samples,p2_samples);
    else
        if sm=='m'
            min_avg_reproj_err=compute_reprojection_error(chosen_phi,epi1_sph,epi2_sph,p1_samples,p2_samples);
        end
    end
end