function [chosen_R,chosen_T,chosen_phi_avg,chosen_phi]=epipole_corrs_to_RT(epi12_sph,p1,p2,ptsnum,varargin)
% Given a pair of epipoles in image sphere and ptsnum
% correspondences in image plane, calculate R,T and phi.
    if ~isempty(varargin)
        K=varargin{1};
        p1=K\p1;
        p2=K\p2;
    end
    p1_samples=p1(:,1:ptsnum);
    p2_samples=p2(:,1:ptsnum);
    epi1_sph=epi12_sph(1:2)';
    epi2_sph=epi12_sph(3:4)';
    epi1_plane=[tan(epi12_sph(1))*cos(epi12_sph(2));tan(epi12_sph(1))*sin(epi12_sph(2));1];
    epi2_plane=[tan(epi12_sph(3))*cos(epi12_sph(4));tan(epi12_sph(3))*sin(epi12_sph(4));1];
    [phi1,phi2,phi1_alt,phi2_alt]=calc_phi_new(p1_samples,p2_samples,epi1_plane,epi2_plane,ptsnum);
    avg_phi1=sum(phi1)/ptsnum;
    avg_phi2=avg_phi1+pi;
    avg_phi1_alt=sum(phi1_alt)/ptsnum;
    avg_phi2_alt=avg_phi1_alt+pi;
    u_gt=p1_samples-repmat(epi1_plane,1,ptsnum);
    ubar_gt=p2_samples-repmat(epi2_plane,1,ptsnum);
    reproj_errs=-ones(8,ptsnum);
    %% Construct R2 and R2_alt.
    v1=epi1_plane/norm(epi1_plane);
    v2=epi2_plane/norm(epi2_plane);
    u2=cross(v1,v2);
    rot1=v1'*v2;
    u2_cross=[0,-u2(3),u2(2);u2(3),0,-u2(1);-u2(2),u2(1),0];
    R2=eye(3)+u2_cross+u2_cross*u2_cross/(1+rot1);
    R2_additional_axis_y=-(epi2_plane(1)+1)/epi2_plane(2);
    R2_additional_cross=[0,-1,R2_additional_axis_y;1,0,-1;-R2_additional_axis_y,1,0];
    R2_additional_cross2=R2_additional_cross*R2_additional_cross;
    R2_additional=eye(3)+2*epi2_plane(2)^2/(2*epi2_plane(2)^2+epi2_plane(1)^2+2*epi2_plane(1)+1)*R2_additional_cross2;
    R2_alt=R2_additional*R2;
    %% Construct R3.
    v2_cross=[0,-v2(3),v2(2);v2(3),0,-v2(1);-v2(2),v2(1),0];
    v2_cross2=v2_cross*v2_cross;
    R3_1=eye(3)+sin(avg_phi1)*v2_cross+(1-cos(avg_phi1))*v2_cross2;
    R3_2=eye(3)+sin(avg_phi2)*v2_cross+(1-cos(avg_phi2))*v2_cross2;
    R3_1_alt=eye(3)+sin(avg_phi1_alt)*v2_cross+(1-cos(avg_phi1_alt))*v2_cross2;
    R3_2_alt=eye(3)+sin(avg_phi2_alt)*v2_cross+(1-cos(avg_phi2_alt))*v2_cross2;
    %% Construct Rs.
    calcR1=R3_1*R2;
    calcR2=R3_2*R2;
    calcR3=R3_1_alt*R2_alt;
    calcR4=R3_2_alt*R2_alt;
    %% Construct Ts.
    T=v2;
    T_alt=-T;
    %% Construct Es.
    T_cross=cross_matrix(T);
    T_alt_cross=cross_matrix(T_alt);
    E1=T_cross*calcR1;
    E2=T_cross*calcR2;
    E3=T_cross*calcR3;
    E4=T_cross*calcR4;
    for m=1:ptsnum
        tan_theta_gt=u_gt(2,m)/u_gt(1,m);
        tan_thetabar1=-(E1(1,2)*tan_theta_gt+E1(1,1))/(E1(2,2)*tan_theta_gt+E1(2,1));
        tan_thetabar2=-(E2(1,2)*tan_theta_gt+E2(1,1))/(E2(2,2)*tan_theta_gt+E2(2,1));
        tan_thetabar3=-(E3(1,2)*tan_theta_gt+E3(1,1))/(E3(2,2)*tan_theta_gt+E3(2,1));
        tan_thetabar4=-(E4(1,2)*tan_theta_gt+E4(1,1))/(E4(2,2)*tan_theta_gt+E4(2,1));
        b1=epi2_plane(2)-tan_thetabar1*epi2_plane(1);
        b2=epi2_plane(2)-tan_thetabar2*epi2_plane(1);
        b3=epi2_plane(2)-tan_thetabar3*epi2_plane(1);
        b4=epi2_plane(2)-tan_thetabar4*epi2_plane(1);
        b5=epi2_plane(2)+tan_thetabar1*epi2_plane(1);
        b6=epi2_plane(2)+tan_thetabar2*epi2_plane(1);
        b7=epi2_plane(2)+tan_thetabar3*epi2_plane(1);
        b8=epi2_plane(2)+tan_thetabar4*epi2_plane(1);
        d1=point_to_line(p2(:,m),epi2_plane,[0;b1;1]);
        d2=point_to_line(p2(:,m),epi2_plane,[0;b2;1]);
        d3=point_to_line(p2(:,m),epi2_plane,[0;b3;1]);
        d4=point_to_line(p2(:,m),epi2_plane,[0;b4;1]);
        d5=point_to_line(p2(:,m),epi2_plane,[0;b5;1]);
        d6=point_to_line(p2(:,m),epi2_plane,[0;b6;1]);
        d7=point_to_line(p2(:,m),epi2_plane,[0;b7;1]);
        d8=point_to_line(p2(:,m),epi2_plane,[0;b8;1]);
        reproj_errs(1,m)=d1;
        reproj_errs(2,m)=d2;
        reproj_errs(3,m)=d3;
        reproj_errs(4,m)=d4;
        reproj_errs(5,m)=d5;
        reproj_errs(6,m)=d6;
        reproj_errs(7,m)=d7;
        reproj_errs(8,m)=d8;
    end
    [~,min_reproj_errs_ind]=min(sum(reproj_errs,2));
    check_depth1=zeros(1,ptsnum);
    check_depth2=zeros(1,ptsnum);
    switch min_reproj_errs_ind
        case 1
            candi_R1=calcR1;
            candi_R2=calcR2;
            for m=1:ptsnum
                check_depth1(m)=v2'*cross(ubar_gt(:,m),cross(candi_R1*u_gt(:,m),v2));
                check_depth2(m)=v2'*cross(ubar_gt(:,m),cross(candi_R2*u_gt(:,m),v2));
            end
            posdepth1=sum(check_depth1<0);
            posdepth2=sum(check_depth2<0);
            switch max([posdepth1,posdepth2])
                    case posdepth1
                        chosen_R=candi_R1;
                        chosen_phi_avg=avg_phi1;
                        chosen_phi=phi1;
                    case posdepth2
                        chosen_R=candi_R2;
                        chosen_phi_avg=avg_phi2;
                        chosen_phi=phi2;
            end
            chosen_T=T;
        case 2
            candi_R1=calcR1;
            candi_R2=calcR2;
            for m=1:ptsnum
                check_depth1(m)=v2'*cross(ubar_gt(:,m),cross(candi_R1*u_gt(:,m),v2));
                check_depth2(m)=v2'*cross(ubar_gt(:,m),cross(candi_R2*u_gt(:,m),v2));
            end
            posdepth1=sum(check_depth1<0);
            posdepth2=sum(check_depth2<0);
            switch max([posdepth1,posdepth2])
                    case posdepth1
                        chosen_R=candi_R1;
                        chosen_phi_avg=avg_phi1;
                        chosen_phi=phi1;
                    case posdepth2
                        chosen_R=candi_R2;
                        chosen_phi_avg=avg_phi2;
                        chosen_phi=phi2;
            end
            chosen_T=T;
        case 3
            candi_R1=calcR3;
            candi_R2=calcR4;
            for m=1:ptsnum
                check_depth1(m)=v2'*cross(ubar_gt(:,m),cross(candi_R1*u_gt(:,m),v2));
                check_depth2(m)=v2'*cross(ubar_gt(:,m),cross(candi_R2*u_gt(:,m),v2));
            end
            posdepth1=sum(check_depth1<0);
            posdepth2=sum(check_depth2<0);
            switch max([posdepth1,posdepth2])
                    case posdepth1
                        chosen_R=candi_R1;
                        chosen_phi_avg=avg_phi1_alt;
                        chosen_phi=phi1_alt;
                    case posdepth2
                        chosen_R=candi_R2;
                        chosen_phi_avg=avg_phi2_alt;
                        chosen_phi=phi2_alt;
            end
            chosen_T=T;
        case 4
            candi_R1=calcR3;
            candi_R2=calcR4;
            for m=1:ptsnum
                check_depth1(m)=v2'*cross(ubar_gt(:,m),cross(candi_R1*u_gt(:,m),v2));
                check_depth2(m)=v2'*cross(ubar_gt(:,m),cross(candi_R2*u_gt(:,m),v2));
            end
            posdepth1=sum(check_depth1<0);
            posdepth2=sum(check_depth2<0);
            switch max([posdepth1,posdepth2])
                    case posdepth1
                        chosen_R=candi_R1;
                        chosen_phi_avg=avg_phi1_alt;
                        chosen_phi=phi1_alt;
                    case posdepth2
                        chosen_R=candi_R2;
                        chosen_phi_avg=avg_phi2_alt;
                        chosen_phi=phi2_alt;
            end
            chosen_T=T;
        case 5
            candi_R1=calcR1;
            candi_R2=calcR2;
            for m=1:ptsnum
                check_depth1(m)=v2'*cross(ubar_gt(:,m),cross(candi_R1*u_gt(:,m),v2));
                check_depth2(m)=v2'*cross(ubar_gt(:,m),cross(candi_R2*u_gt(:,m),v2));
            end
            posdepth1=sum(check_depth1<0);
            posdepth2=sum(check_depth2<0);
            switch max([posdepth1,posdepth2])
                    case posdepth1
                        chosen_R=candi_R1;
                        chosen_phi_avg=avg_phi1;
                        chosen_phi=phi1;
                    case posdepth2
                        chosen_R=candi_R2;
                        chosen_phi_avg=avg_phi2;
                        chosen_phi=phi2;
            end
            chosen_T=T_alt;
        case 6
            candi_R1=calcR1;
            candi_R2=calcR2;
            for m=1:ptsnum
                check_depth1(m)=v2'*cross(ubar_gt(:,m),cross(candi_R1*u_gt(:,m),v2));
                check_depth2(m)=v2'*cross(ubar_gt(:,m),cross(candi_R2*u_gt(:,m),v2));
            end
            posdepth1=sum(check_depth1<0);
            posdepth2=sum(check_depth2<0);
            switch max([posdepth1,posdepth2])
                    case posdepth1
                        chosen_R=candi_R1;
                        chosen_phi_avg=avg_phi1;
                        chosen_phi=phi1;
                    case posdepth2
                        chosen_R=candi_R2;
                        chosen_phi_avg=avg_phi2;
                        chosen_phi=phi2;
            end
            chosen_T=T_alt;
        case 7
            candi_R1=calcR3;
            candi_R2=calcR4;
            for m=1:ptsnum
                check_depth1(m)=v2'*cross(ubar_gt(:,m),cross(candi_R1*u_gt(:,m),v2));
                check_depth2(m)=v2'*cross(ubar_gt(:,m),cross(candi_R2*u_gt(:,m),v2));
            end
            posdepth1=sum(check_depth1<0);
            posdepth2=sum(check_depth2<0);
            switch max([posdepth1,posdepth2])
                    case posdepth1
                        chosen_R=candi_R1;
                        chosen_phi_avg=avg_phi1_alt;
                        chosen_phi=phi1_alt;
                    case posdepth2
                        chosen_R=candi_R2;
                        chosen_phi_avg=avg_phi2_alt;
                        chosen_phi=phi2_alt;
            end
            chosen_T=T_alt;
        case 8
            candi_R1=calcR3;
            candi_R2=calcR4;
            for m=1:ptsnum
                check_depth1(m)=v2'*cross(ubar_gt(:,m),cross(candi_R1*u_gt(:,m),v2));
                check_depth2(m)=v2'*cross(ubar_gt(:,m),cross(candi_R2*u_gt(:,m),v2));
            end
            posdepth1=sum(check_depth1<0);
            posdepth2=sum(check_depth2<0);
            switch max([posdepth1,posdepth2])
                    case posdepth1
                        chosen_R=candi_R1;
                        chosen_phi_avg=avg_phi1_alt;
                        chosen_phi=phi1_alt;
                    case posdepth2
                        chosen_R=candi_R2;
                        chosen_phi_avg=avg_phi2_alt;
                        chosen_phi=phi2_alt;
            end
            chosen_T=T_alt;
    end
end