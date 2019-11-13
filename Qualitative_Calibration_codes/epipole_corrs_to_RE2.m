function RE=epipole_corrs_to_RE2(epi12_sph,p1,p2,ptsnum,sm,varargin)
    if ~isempty(varargin)
        K1=varargin{1};
        K2=varargin{2};
        [R,T]=epipole_corrs_to_RT2(epi12_sph,p1,p2,ptsnum,K1,K2);
    else
        [R,T]=epipole_corrs_to_RT2(epi12_sph,p1,p2,ptsnum);
    end
    p1_samples=p1(:,1:ptsnum);
    p2_samples=p2(:,1:ptsnum);
%     reproj_errs=RT_to_RE(R,T,p1_samples,p2_samples,ptsnum);
    if isempty(varargin{1})
        reproj_errs=compute_reprojection_error_all(R,T,p1_samples,p2_samples,ptsnum);
    else
        K1=varargin{1};
        K2=varargin{2};
        reproj_errs=compute_reprojection_error_all(R,T,p1_samples,p2_samples,ptsnum,K1,K2);
    end
    if sm=='s'
        RE=sum(reproj_errs)/ptsnum;
    else
        if sm=='m'
            RE=reproj_errs;
        end
    end
end