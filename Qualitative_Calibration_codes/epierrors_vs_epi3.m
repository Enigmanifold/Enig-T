trialnum=500;
total_total_cell_sph=cell(trialnum,1);
for a=1:trialnum
    run('generate_RT.m');
    e1=R'*T;
    e1=e1/e1(3);
    e2=T/T(3);
    Theta_seg=11;
    Phi_seg=10;
    Theta_lim=pi/2*(1-1/Theta_seg);
    Theta=linspace(0,Theta_lim,Theta_seg);
    Phi=linspace(0,Theta_lim,Theta_seg);
    dTheta=Theta(2);
    dPhi=Phi(2);
    e1Theta=acos(1/norm(e1));
    e2Theta=acos(1/norm(e2));
    e1Phi=atan(e1(2)/e1(1));
    if e1(1)<0
        e1Phi=e1Phi+pi;
    end
    e2Phi=atan(e2(2)/e2(1));
    if e2(1)<0
        e2Phi=e2Phi+pi;
    end
    e1_sph=[e1Theta;e1Phi];
    e2_sph=[e2Theta;e2Phi];
    ptsnum=size(DPoints,2);
    ptsnum=30;
    Thetas=[0,-dTheta,dTheta];
    Phis=[0,-dPhi,dPhi];
    e1s_sph=[Thetas;Phis]+repmat(e1_sph,1,3);
    e2s_sph=[Thetas;Phis]+repmat(e2_sph,1,3);
    epi1s_sph=ones(2,9);
    epi2s_sph=ones(2,9);
    epi1s_sph(2,:)=repmat(e1s_sph(2,:),1,3);
    epi2s_sph(2,:)=repmat(e2s_sph(2,:),1,3);
    total_cell_sph=cell(81,8);
    counter=1;
    for m=1:3
        for n=1:3
            epi1s_sph(1,3*m-2:3*m)=e1s_sph(1,m);
            epi2s_sph(1,3*m-2:3*m)=e2s_sph(1,m);
        end
    end
    for m=1:9
        for n=1:9
            epi1_sph=epi1s_sph(:,m);
            epi2_sph=epi2s_sph(:,n);
            if epi1_sph(1)>=pi/2 || epi2_sph(1)>=pi/2
                continue
            end
            epi1=[sin(epi1_sph(1))*cos(epi1_sph(2))*sec(epi1_sph(1));sin(epi1_sph(1))*sin(epi1_sph(2))*sec(epi1_sph(1));1];
            epi2=[sin(epi2_sph(1))*cos(epi2_sph(2))*sec(epi2_sph(1));sin(epi2_sph(1))*sin(epi2_sph(2))*sec(epi2_sph(1));1];
            [phi1,phi2,phi1_alt,phi2_alt]=calc_phi(p1,p2,epi1,epi2,ptsnum);
            p1_to_e1=p1(:,1:ptsnum)-repmat(epi1,1,ptsnum);
            p2_to_e2=p2(:,1:ptsnum)-repmat(epi2,1,ptsnum);
            u0=p1_to_e1./vecnorm(p1_to_e1);
            u0_bar=p2_to_e2./vecnorm(p2_to_e2);
            v1=epi1./norm(epi1);
            v2=epi2./norm(epi2);
            v2_cross=[0,-v2(3),v2(2);v2(3),0,-v2(1);-v2(2),v2(1),0];
            v2_cross2=v2_cross*v2_cross;
            v1_cross=[0,-v1(3),v1(2);v1(3),0,-v1(1);-v1(2),v1(1),0];
            v1_cross2=v1_cross*v1_cross;
            u1=cross(v1,v2);
            rot1=v1'*v2;
            u1_cross=[0,-u1(3),u1(2);u1(3),0,-u1(1);-u1(2),u1(1),0];
            R1=eye(3)+u1_cross+u1_cross*u1_cross/(1+rot1);
            R1_ambi_y=-(e2(1)+1)/e2(2);
            R1_ambi_axis=[1;R1_ambi_y;1]/norm([1,R1_ambi_y,1]);
            R1_ambi_axis_cross=[0,-R1_ambi_axis(3),R1_ambi_axis(2);R1_ambi_axis(3),0,-R1_ambi_axis(1);-R1_ambi_axis(2),R1_ambi_axis(1),0];
            R1_ambi=eye(3)+2*(R1_ambi_axis_cross*R1_ambi_axis_cross);
            R1_alt=R1_ambi*R1;
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
            for o=1:ptsnum
                R2_1{o}=eye(3)+sinphi1(o)*v2_cross+(1-cosphi1(o))*v2_cross2;
                R2_2{o}=eye(3)+sinphi2(o)*v2_cross+(1-cosphi2(o))*v2_cross2;
                R2_3{o}=eye(3)+sinphi1_alt(o)*v2_cross+(1-cosphi1_alt(o))*v2_cross2;
                R2_4{o}=eye(3)+sinphi2_alt(o)*v2_cross+(1-cosphi2_alt(o))*v2_cross2;
            end
            calcR1_all=cell(1,ptsnum);
            calcR2_all=cell(1,ptsnum);
            calcR3_all=cell(1,ptsnum);
            calcR4_all=cell(1,ptsnum);
            for o=1:ptsnum
                calcR1_all{o}=R2_1{o}*R1;
                calcR2_all{o}=R2_2{o}*R1;
                calcR3_all{o}=R2_3{o}*R1_alt;
                calcR4_all{o}=R2_4{o}*R1_alt;
            end
            calcR1234_all={calcR1_all,calcR2_all,calcR3_all,calcR4_all};
            calcR_all=cell(1,ptsnum);
            for o=1:4
                axs=zeros(3,ptsnum);
                angs=zeros(1,ptsnum);
                calcR_tmp=calcR1234_all{o};
                for p=1:ptsnum
                    axang=rotm2axang(calcR_tmp{p});
                    axs(:,p)=axang(1:3)';
                    angs(p)=axang(4);
                end
                ax=mean(axs,2);
                ax=ax/norm(ax);
                ang=mean(angs);
                ax_cross=[0,-ax(3),ax(2);ax(3),0,-ax(1);-ax(2),ax(1),0];
                ax_cross2=ax_cross*ax_cross;
                calcR_all{o}=eye(3)+sin(ang)*ax_cross+(1-cos(ang))*ax_cross2;
            end
            check_depth1=zeros(1,ptsnum);
            check_depth2=zeros(1,ptsnum);
            check_depth3=zeros(1,ptsnum);
            check_depth4=zeros(1,ptsnum);
            for o=1:ptsnum
                check_depth1(o)=v2'*cross(u0_bar(:,o),cross(calcR_all{1}*u0(:,o),v2));
                check_depth2(o)=v2'*cross(u0_bar(:,o),cross(calcR_all{2}*u0(:,o),v2));
                check_depth3(o)=v2'*cross(u0_bar(:,o),cross(calcR_all{3}*u0(:,o),v2));
                check_depth4(o)=v2'*cross(u0_bar(:,o),cross(calcR_all{4}*u0(:,o),v2));
            end
            posdepth1=sum(check_depth1<0);
            posdepth2=sum(check_depth2<0);
            posdepth3=sum(check_depth3<0);
            posdepth4=sum(check_depth4<0);
            mpdinds=find([posdepth1,posdepth2,posdepth3,posdepth4]==max([posdepth1,posdepth2,posdepth3,posdepth4]));
            if length(mpdinds)==1
                switch mpdinds
                    case 1
                        calcR=calcR_all{1};
                        ind2=1;
                        phis2=phi1;
                    case 2
                        calcR=calcR_all{2};
                        ind2=2;
                        phis2=phi2;
                    case 3
                        calcR=calcR_all{3};
                        ind2=3;
                        phis2=phi1_alt;
                    case 4
                        calcR=calcR_all{4};
                        ind2=4;
                        phis2=phi2_alt;
                end
            else
                phis=[phi1;phi2;phi1_alt;phi2_alt];
                phis2_std=-ones(length(mpdinds),1);
                for o=1:length(mpdinds)                   
                    phis2_candi=phis(mpdinds(o),:);
                    phis2_std(o)=std(phis2_candi);
                end
                [~,phis2_std_min]=min(phis2_std);
                phis2=phis(mpdinds(phis2_std_min),:);
            end
            phis2_c0=-ones(1,ptsnum);
            for o=1:ptsnum
                [~,phi2_c0_ind]=min([abs(phis2(o)),abs(phis2(o)-pi),abs(phis2(o)-2*pi)]);
                switch phi2_c0_ind
                    case 1
                        phis2_c0(o)=phis2(o);
                    case 2
                        phis2_c0(o)=phis2(o)-pi;
                    case 3
                        phis2_c0(o)=phis2(o)-2*pi;
                end
            end
            minphistd2=std(phis2_c0.*180./pi);
            total_cell_sph{counter,1}=epi1_sph'.*180./pi;
            total_cell_sph{counter,2}=epi2_sph'.*180./pi;
            total_cell_sph{counter,3}=(epi1_sph'-e1_sph').*180./pi;
            total_cell_sph{counter,4}=(epi2_sph'-e2_sph').*180./pi;
            total_cell_sph{counter,5}=epi1';
            total_cell_sph{counter,6}=epi2';
    %         total_cell_sph{counter,3}=phis;
            total_cell_sph{counter,7}=phis2.*180/pi;
    %         total_cell_sph{counter,4}=minphistd;
            total_cell_sph{counter,8}=minphistd2;
            total_cell_sph{counter,9}=T_scale;
            total_cell_sph{counter,10}=[posdepth1,posdepth2,posdepth3,posdepth4];
            counter=counter+1;
        end
    end
    total_total_cell_sph{a}=total_cell_sph;
end