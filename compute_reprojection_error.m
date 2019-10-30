function min_reproj_errs=compute_reprojection_error(phi,epi1_sph,epi2_sph,p1,p2)
    ptsnum=size(p1,2);
    reproj_errs=-ones(4,ptsnum);
    [epi1_plane(1),epi1_plane(2)]=imsphere2implane(epi1_sph(1),epi1_sph(2));
    [epi2_plane(1),epi2_plane(2)]=imsphere2implane(epi2_sph(1),epi2_sph(2));
    epi1_plane=[epi1_plane';1];
    epi2_plane=[epi2_plane';1];
    [candi_R,candi_R_alt,candi_T,candi_T_alt]=epipoles_phi_to_RT(epi1_plane,epi2_plane,phi,1);
    candi_T_cross=cross_matrix(candi_T);
    candi_T_alt_cross=cross_matrix(candi_T_alt);
    E1=candi_T_cross*candi_R;
    E2=candi_T_cross*candi_R_alt;
    E3=candi_T_alt_cross*candi_R;
    E4=candi_T_alt_cross*candi_R_alt;
    u=p1-repmat(epi1_plane,1,ptsnum);
    for m=1:ptsnum
        tan_theta=u(2,m)/u(1,m);
        tan_thetabar1=-(E1(1,2)*tan_theta+E1(1,1))/(E1(2,2)*tan_theta+E1(2,1));
        tan_thetabar2=-(E2(1,2)*tan_theta+E2(1,1))/(E2(2,2)*tan_theta+E2(2,1));
        tan_thetabar3=-(E3(1,2)*tan_theta+E3(1,1))/(E3(2,2)*tan_theta+E3(2,1));
        tan_thetabar4=-(E4(1,2)*tan_theta+E4(1,1))/(E4(2,2)*tan_theta+E4(2,1));
        b1=epi2_plane(2)-tan_thetabar1*epi2_plane(1);
        b2=epi2_plane(2)-tan_thetabar2*epi2_plane(1);
        b3=epi2_plane(2)-tan_thetabar3*epi2_plane(1);
        b4=epi2_plane(2)-tan_thetabar4*epi2_plane(1);
        d1=point_to_line(p2(:,m),epi2_plane,[0;b1;1]);
        d2=point_to_line(p2(:,m),epi2_plane,[0;b2;1]);
        d3=point_to_line(p2(:,m),epi2_plane,[0;b3;1]);
        d4=point_to_line(p2(:,m),epi2_plane,[0;b4;1]);
        reproj_errs(1,m)=d1;
        reproj_errs(2,m)=d2;
        reproj_errs(3,m)=d3;
        reproj_errs(4,m)=d4;
    end
    [~,min_reproj_errs_ind]=min(sum(reproj_errs,2));
    min_reproj_errs=reproj_errs(min_reproj_errs_ind,:);
end
        