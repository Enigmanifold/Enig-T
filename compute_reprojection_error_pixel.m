function errs=compute_reprojection_error_pixel(R,T,K,p1_pixel,p2_pixel,ptsnum)
    p1_pixel=p1_pixel(:,1:ptsnum);
    p2_pixel=p2_pixel(:,1:ptsnum);
    T_cross=cross_matrix(T);
    E=T_cross*R;
    invK=inv(K);
    errs = zeros(1,ptsnum);
    F = invK'*E*invK;
    abc = F*p1_pixel;
    ab = abc(1:2,:);
    norm_const = vecnorm(ab,2,1);
    p2_pixel_trans=p2_pixel';
    for n = 1:ptsnum
        err = abs(p2_pixel_trans(n,:)*abc(:,n));
        errs(1,n) = err;
    end
    errs = errs./norm_const;
end