function errs=epiphi_to_RE(epi12_sph_phi,p1,p2,ptsnum,varargin)
% Find the pair of epipoles and phi that gives the smallest reprojection
% error.
    narginchk(4,5);
    if nargin==5
        K=varargin{1};
    end
    p1_samples=p1(:,1:ptsnum);
    p2_samples=p2(:,1:ptsnum);
    epi1_plane=[tan(epi12_sph_phi(1))*cos(epi12_sph_phi(2));tan(epi12_sph_phi(1))*sin(epi12_sph_phi(2));1];
    epi2_plane=[tan(epi12_sph_phi(3))*cos(epi12_sph_phi(4));tan(epi12_sph_phi(3))*sin(epi12_sph_phi(4));1];
    phi=epi12_sph_phi(5);
    [R,R_alt,T,T_alt]=epipoles_phi_to_RT(epi1_plane,epi2_plane,phi,-1);
    if isempty(varargin)
%         disp('Focal length regime.')
        errs1=compute_reprojection_error_all(R,T,p1_samples,p2_samples,ptsnum);
        errs2=compute_reprojection_error_all(R_alt,T,p1_samples,p2_samples,ptsnum);
        errs3=compute_reprojection_error_all(R,T_alt,p1_samples,p2_samples,ptsnum);
        errs4=compute_reprojection_error_all(R_alt,T_alt,p1_samples,p2_samples,ptsnum);
    else
%         disp('Pixel regime.')
        errs1=compute_reprojection_error_all(R,T,p1_samples,p2_samples,ptsnum,K);
        errs2=compute_reprojection_error_all(R_alt,T,p1_samples,p2_samples,ptsnum,K);
        errs3=compute_reprojection_error_all(R,T_alt,p1_samples,p2_samples,ptsnum,K);
        errs4=compute_reprojection_error_all(R_alt,T_alt,p1_samples,p2_samples,ptsnum,K);
    end
    errs_all=[errs1;errs2;errs3;errs4];
    errs1_norm=norm(errs1);
    errs2_norm=norm(errs2);
    errs3_norm=norm(errs3);
    errs4_norm=norm(errs4);
    [~,errs_min_ind]=min([errs1_norm,errs2_norm,errs3_norm,errs4_norm]);
    errs=errs_all(errs_min_ind,:);
end