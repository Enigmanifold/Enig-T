function errs=compute_reprojection_error_fl(R,T,p1,p2,ptsnum)
    p1=p1(:,1:ptsnum);
    p2=p2(:,1:ptsnum);
    T_cross=cross_matrix(T);
    E=T_cross*R;
    errs = zeros(1,ptsnum);
    abc = E*p1;
    ab = abc(1:2,:);
    norm_const = vecnorm(ab,2,1);
    p2_trans=p2';
    for n = 1:ptsnum
        err = abs(p2_trans(n,:)*abc(:,n));
        errs(1,n) = err;
    end
    errs = errs./norm_const;
end