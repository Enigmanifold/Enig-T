%% Getting epipoles and theta.
run('generate_RT.m');
Rinv=R';
e1=R'*T;
e1=e1/e1(3);
e2=T/T(3);
ptsnum=size(DPoints,2);
ptsnum=20;
radii=[0,0.1,0.3,1,2,3,5,7,10,20,30,45,60];
angles=[0,pi/6,2*pi/6,3*pi/6,4*pi/6,5*pi/6,pi,7*pi/6,8*pi/6,9*pi/6,10*pi/6,11*pi/6];
axs=zeros(3,ptsnum);
angs=zeros(1,ptsnum);
e1_epis=ones(3,size(radii,2),size(angles,2));
e2_epis=ones(3,size(radii,2),size(angles,2));
counter=1;
for m=1:size(radii,2)
    for n=1:size(angles,2)
        e1_epis(:,m,n)=[e1(1)+radii(m)*cos(angles(n));e1(2)+radii(m)*sin(angles(n));1];
        e2_epis(:,m,n)=[e2(1)+radii(m)*cos(angles(n));e2(2)+radii(m)*sin(angles(n));1];
    end
end
errormap=-ones(size(radii,2),size(angles,2),size(radii,2),size(angles,2));
for a=1:size(errormap,1)
    for b=1:size(errormap,2)
        for c=1:size(errormap,3)
            for d=1:size(errormap,4)
                p1_to_e1=p1(:,1:ptsnum)-repmat(e1_epis(:,a,b),1,ptsnum);
                p2_to_e2=p2(:,1:ptsnum)-repmat(e2_epis(:,c,d),1,ptsnum);
                u0=p1_to_e1./vecnorm(p1_to_e1);
                u0_bar=p2_to_e2./vecnorm(p2_to_e2);
                theta1=zeros(1,ptsnum);
                theta2=zeros(1,ptsnum);
                for m=1:ptsnum
                    th1=atan(p1_to_e1(2,m)/p1_to_e1(1,m));
                    th2=atan(p2_to_e2(2,m)/p2_to_e2(1,m));
                    if p1_to_e1(1,m)<0
                        th1=th1+pi;
                    end
                    if th1<0
                        th1=th1+2*pi;
                    end
                    if p2_to_e2(1,m)<0
                        th2=th2+pi;
                    end
                    if th2<0
                        th2=th2+2*pi;
                    end
                    theta1(m)=th1;
                    theta2(m)=th2;
                end
                %% Calculate R1(aligning two epipoles).
                v1=e1/norm(e1);
                v2=e2/norm(e2);
                u1=cross(v1,v2);
                % u1=u1/norm(u1);
                rot1=v1'*v2;
                u1_cross=[0,-u1(3),u1(2);u1(3),0,-u1(1);-u1(2),u1(1),0];
                u1_cross_alt=-u1_cross;
                % R1=eye(3)+sin(rot1)*u1_cross+(1-cos(rot1))*(u1_cross*u1_cross);
                R1=eye(3)+u1_cross+u1_cross*u1_cross/(1+rot1);
                R1_ambi_y=-(e2(1)+1)/e2(2);
                R1_ambi_axis=[1;R1_ambi_y;1]/norm([1,R1_ambi_y,1]);
                R1_ambi_axis_cross=[0,-R1_ambi_axis(3),R1_ambi_axis(2);R1_ambi_axis(3),0,-R1_ambi_axis(1);-R1_ambi_axis(2),R1_ambi_axis(1),0];
                R1_ambi=eye(3)+2*(R1_ambi_axis_cross*R1_ambi_axis_cross);
                R1_alt=R1_ambi*R1;
                %% Calculate R2(making epipolar half-lines coplanar).
                N1_array=zeros(3,ptsnum);
                N2_array=zeros(3,ptsnum);
                N2_array_alt=zeros(3,ptsnum);
                for m=1:ptsnum
                    N1=cross(u0(:,m),v1);
                    N1=N1/norm(N1);
                    N2=R1*N1;
                    N2_alt=R1_ambi*N2;
                    N1_array(:,m)=N1;
                    N2_array(:,m)=N2;
                    N2_array_alt(:,m)=N2_alt;
                end

                v2_cross=[0,-v2(3),v2(2);v2(3),0,-v2(1);-v2(2),v2(1),0];
                v2_cross2=v2_cross*v2_cross;
                v1_cross=[0,-v1(3),v1(2);v1(3),0,-v1(1);-v1(2),v1(1),0];
                v1_cross2=v1_cross*v1_cross;
                vn=zeros(3,ptsnum);
                un=zeros(1,ptsnum);
                uvn=zeros(1,ptsnum);
                vn_alt=zeros(3,ptsnum);
                un_alt=zeros(1,ptsnum);
                uvn_alt=zeros(1,ptsnum);
                phi1=zeros(1,ptsnum);
                phi2=zeros(1,ptsnum);
                phi1_alt=zeros(1,ptsnum);
                phi2_alt=zeros(1,ptsnum);
                for m=1:ptsnum
                    vn=v2_cross*N2_array(:,m);
                    vn_alt=v2_cross*N2_array_alt(:,m);
                    %     v2n(:,m)=v2_cross2*N2_array(:,m);
                    %     v2n_alt(:,m)=v2_cross2*N2_array_alt(:,m);
                    un=u0_bar(:,m)'*N2_array(:,m);
                    un_alt=u0_bar(:,m)'*N2_array_alt(:,m);
                    uvn=u0_bar(:,m)'*vn;
                    uvn_alt=u0_bar(:,m)'*vn_alt;
                    tanphi=-un/uvn;
                    tanphi_alt=-un_alt/uvn_alt;
                    ph1=atan(tanphi);
                    ph1_alt=atan(tanphi_alt);
                    if ph1<0
                        ph1=ph1+pi;
                    end
                    if ph1_alt<0
                        ph1_alt=ph1_alt+pi;
                    end
                    ph2=ph1+pi;
                    ph2_alt=ph1_alt+pi;
                    phi1(m)=ph1;
                    phi2(m)=ph2;
                    phi1_alt(m)=ph1_alt;
                    phi2_alt(m)=ph2_alt;
                end
                sinphi1=sin(phi1);
                cosphi1=cos(phi1);
                sinphi2=sin(phi2);
                cosphi2=cos(phi2);
                sinphi1_alt=sin(phi1_alt);
                cosphi1_alt=cos(phi1_alt);
                sinphi2_alt=sin(phi2_alt);
                cosphi2_alt=cos(phi2_alt);
                R2_1=cell(1,ptsnum);
                R2_2=cell(1,ptsnum);
                R2_3=cell(1,ptsnum);
                R2_4=cell(1,ptsnum);
                for m=1:ptsnum
                    R2_1{m}=eye(3)+sinphi1(m)*v2_cross+(1-cosphi1(m))*v2_cross2;
                    R2_2{m}=eye(3)+sinphi2(m)*v2_cross+(1-cosphi2(m))*v2_cross2;
                    R2_3{m}=eye(3)+sinphi1_alt(m)*v2_cross+(1-cosphi1_alt(m))*v2_cross2;
                    R2_4{m}=eye(3)+sinphi2_alt(m)*v2_cross+(1-cosphi2_alt(m))*v2_cross2;
                end
                pi_around_T=eye(3)+2*v2_cross2;
                calcR1_all=cell(1,ptsnum);
                calcR2_all=cell(1,ptsnum);
                calcR3_all=cell(1,ptsnum);
                calcR4_all=cell(1,ptsnum);
                for m=1:ptsnum
                    calcR1_all{m}=R2_1{m}*R1;
                    calcR2_all{m}=R2_2{m}*R1;
                    calcR3_all{m}=R2_3{m}*R1_alt;
                    calcR4_all{m}=R2_4{m}*R1_alt;
                end
%                 phi1_std=std(phi1);
%                 phi2_std=std(phi2);
%                 phi1_alt_std=std(phi1_alt);
%                 phi2_alt_std=std(phi2_alt);
%                 switch min([phi1_std,phi2_std,phi1_alt_std,phi2_alt_std])
%                     case phi1_std
%                         calcR_all=calcR1_all;
%                     case phi2_std
%                         calcR_all=calcR2_all;
%                     case phi1_alt_std
%                         calcR_all=calcR3_all;
%                     case phi2_alt_std
%                         calcR_all=calcR4_all;
%                 end
                calcR1234_all={calcR1_all,calcR2_all,calcR3_all,calcR4_all};
                calcR_all=cell(1,ptsnum);
                for m=1:4
                    axs=zeros(3,ptsnum);
                    angs=zeros(1,ptsnum);
                    calcR_tmp=calcR1234_all{m};
                    for n=1:ptsnum
                        axang=rotm2axang(calcR_tmp{n});
                        axs(:,n)=axang(1:3)';
                        angs(n)=axang(4);
                    end
                    ax=mean(axs,2);
                    ax=ax/norm(ax);
                    ang=mean(angs);
                    ax_cross=[0,-ax(3),ax(2);ax(3),0,-ax(1);-ax(2),ax(1),0];
                    ax_cross2=ax_cross*ax_cross;
                    calcR_all{m}=eye(3)+sin(ang)*ax_cross+(1-cos(ang))*ax_cross2;
                end
                check_depth1=zeros(1,ptsnum);
                check_depth2=zeros(1,ptsnum);
                check_depth3=zeros(1,ptsnum);
                check_depth4=zeros(1,ptsnum);
                for m=1:ptsnum
                    check_depth1(m)=v2'*cross(u0_bar(:,m),cross(calcR_all{1}*u0(:,m),v2));
                    check_depth2(m)=v2'*cross(u0_bar(:,m),cross(calcR_all{2}*u0(:,m),v2));
                    check_depth3(m)=v2'*cross(u0_bar(:,m),cross(calcR_all{3}*u0(:,m),v2));
                    check_depth4(m)=v2'*cross(u0_bar(:,m),cross(calcR_all{4}*u0(:,m),v2));
                end
                posdepth1=sum(check_depth1<0);
                posdepth2=sum(check_depth2<0);
                posdepth3=sum(check_depth3<0);
                posdepth4=sum(check_depth4<0);
                switch max([posdepth1,posdepth2,posdepth3,posdepth4])
                    case posdepth1
                        calcR=calcR_all{1};
                        calcphi=phi1;
                    case posdepth2
                        calcR=calcR_all{2};
                        calcphi=phi2;
                    case posdepth3
                        calcR=calcR_all{3};
                        calcphi=phi1_alt;
                    case posdepth4
                        calcR=calcR_all{4};
                        calcphi=phi2_alt;
                end
                error=acos((trace(Rinv*calcR)-1)/2);
                errormap(a,b,c,d)=error;
            end
        end
    end
end


    %     error=-ones(1,ptsnum);
    %     ind=-ones(1,ptsnum);
    %     Rinv=R';
    %     for m=1:ptsnum
    %         if isreal(calcR1{m})
    %             R1error=acos((trace(Rinv*calcR1{m})-1)/2);
    %         else
    %             R1error=99;
    %         end
    %         if isreal(calcR2{m})
    %             R2error=acos((trace(Rinv*calcR2{m})-1)/2);
    %         else
    %             R2error=99;
    %         end
    %         if isreal(calcR3{m})
    %             R3error=acos((trace(Rinv*calcR3{m})-1)/2);
    %         else
    %             R3error=99;
    %         end
    %         if isreal(calcR4{m})
    %             R4error=acos((trace(Rinv*calcR4{m})-1)/2);
    %         else
    %             R4error=99;
    %         end
    %         [error(m),ind(m)]=min([R1error,R2error,R3error,R4error]);
    %     end
    %     if min(ind)==max(ind)
    %         disp(strcat('From calcR',num2str(min(ind))));
    %         switch min(ind)
    %             case 1
    %                 calcR=calcR1{1};
    %             case 2
    %                 calcR=calcR2{1};
    %             case 3
    %                 calcR=calcR3{1};
    %             case 4
    %                 calcR=calcR4{1};
    %         end
    %     else
    %         disp('Indeterminate R!');
    %     end
    %     % rotate_calcR=pi_around_T*calcR;
    %     % epi1=rotate_calcR'*T;
    %     % epi1=epi1/epi1(3);
    % 
    %     check_depth1=zeros(1,ptsnum);
    %     check_depth2=zeros(1,ptsnum);
    %     check_depth3=zeros(1,ptsnum);
    %     check_depth4=zeros(1,ptsnum);
    %     for m=1:ptsnum
    %         check_depth1(m)=v2'*cross(u0_bar(:,m),cross(calcR1{m}*u0(:,m),v2));
    %         check_depth2(m)=v2'*cross(u0_bar(:,m),cross(calcR2{m}*u0(:,m),v2));
    %         check_depth3(m)=v2'*cross(u0_bar(:,m),cross(calcR3{m}*u0(:,m),v2));
    %         check_depth4(m)=v2'*cross(u0_bar(:,m),cross(calcR4{m}*u0(:,m),v2));
    %     end
    %     disp(strcat('Total number of points: ',num2str(ptsnum)));
    %     disp(strcat('Number of positive depths for calcR1 is ',num2str(sum(check_depth1<0))));
    %     disp(strcat('Number of positive depths for calcR2 is ',num2str(sum(check_depth2<0))));
    %     disp(strcat('Number of positive depths for calcR3 is ',num2str(sum(check_depth3<0))));
    %     disp(strcat('Number of positive depths for calcR4 is ',num2str(sum(check_depth4<0))));
        % theta_diff=theta2-theta1;
        % for m=1:ptsnum
        %     if theta_diff(m)>0
        %         theta_diff(m)=theta_diff(m)-2*pi;
        %     end
        % end
        % plot(linspace(0,ptsnum-1,ptsnum),theta_diff)