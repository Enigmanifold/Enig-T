function min_avg_reproj_err=compute_reprojection_error_avg(phi,epi1_sph,epi2_sph,p1,p2)
    min_reproj_errs=compute_reprojection_error(phi,epi1_sph,epi2_sph,p1,p2);
    min_avg_reproj_err=sum(min_reproj_errs)/length(min_reproj_errs);
end