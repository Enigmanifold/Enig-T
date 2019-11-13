function reproj_errs=RT_to_RE(R,T,p1,p2,ptsnum)
    e1=R'*T;
    e1=e1/e1(3);
    e2=T/T(3);
    p1=p1(:,1:ptsnum);
    p2=p2(:,1:ptsnum);
    T_cross=cross_matrix(T);
    E=T_cross*R;
    u=p1-repmat(e1,1,ptsnum);
    reproj_errs=-ones(1,ptsnum);
    for m=1:ptsnum
        tan_theta=u(2,m)/u(1,m);
        tan_thetabar=-(E(1,2)*tan_theta+E(1,1))/(E(2,2)*tan_theta+E(2,1));
        b1=e2(2)-tan_thetabar*e2(1);
        d1=point_to_line(p2(:,m),e2,[0;b1;1]);
        reproj_errs(1,m)=d1;
    end
end