function RE=epipole_corrs_to_RE2(epi12_sph,p1,p2,ptsnum,sm,varargin)
    [R,T,chosenphi]=epipole_corrs_to_RT(epi12_sph,p1,p2,ptsnum);
    p1_samples=p1(:,1:ptsnum);
    p2_samples=p2(:,1:ptsnum);
%     reproj_errs=RT_to_RE(R,T,p1_samples,p2_samples,ptsnum);
    reproj_errs=compute_reprojection_error_all(R,T,p1_samples,p2_samples,ptsnum);
    if sm=='s'
        RE=sum(reproj_errs)/ptsnum;
    else
        if sm=='m'
            RE=reproj_errs;
        end
    end
end