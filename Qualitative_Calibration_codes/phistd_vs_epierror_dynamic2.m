%% Errors due to increasing constant distances.
run('generate_RT.m');
Rinv=R';
e1=R'*T;
e1=e1/e1(3);
e2=T/T(3);
ptsnum=size(DPoints,2);
ptsnum=20;
% randr=normrnd(0,1,[2,ptsnum]); % Introduce noise.
% randtheta=rand(2,ptsnum).*2.*pi;
% p1noise=[randr(1,:)+cos(randtheta(1,:));randr(1,:)+sin(randtheta(1,:));zeros(1,ptsnum)];
% p2noise=[randr(2,:)+cos(randtheta(2,:));randr(2,:)+sin(randtheta(2,:));zeros(1,ptsnum)];
% p1(:,1:ptsnum)=p1(:,1:ptsnum)+p1noise;
% p2(:,1:ptsnum)=p2(:,1:ptsnum)+p2noise;
mu=zeros(1,2);
t=linspace(-50,50,100).*0.1;
tlen=length(t);
errs=exp(t);
errs(1)=0;
samplen=100;
minstdarray=-ones(tlen,samplen);
minstdmultiarray=-ones(tlen,samplen,tlen);
for f=1:tlen
    rndangles1=rand(1,1).*2.*pi;
    for a=1:tlen
        err=errs(a);
%         rndangles1=rand(1,samplen).*2.*pi;
        rndangles2=linspace(0,samplen-1,samplen)./50.*pi;
%         rndepi1err=[errs(f).*cos(rndangles1);errs(f).*sin(rndangles1);zeros(1,samplen)];
        rndepi1err=repmat([errs(f).*cos(rndangles1);errs(f).*sin(rndangles1);0],1,samplen);
        rndepi2err=[err.*cos(rndangles2);err.*sin(rndangles2);zeros(1,samplen)];
        epi1_samples=repmat(e1,1,samplen)+rndepi1err;
        epi2_samples=repmat(e2,1,samplen)+rndepi2err;
        for b=1:samplen
            p1_to_e1=p1(:,1:ptsnum)-repmat(epi1_samples(:,b),1,ptsnum);
            p2_to_e2=p2(:,1:ptsnum)-repmat(epi2_samples(:,b),1,ptsnum);
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
            phi1_std=std(phi1);
            phi2_std=std(phi2);
            phi1_alt_std=std(phi1_alt);
            phi2_alt_std=std(phi2_alt);
            minstd=min([phi1_std,phi2_std,phi1_alt_std,phi2_alt_std]);
            minstdarray(a,b)=minstd;
        end
    end
    minstdmultiarray(:,:,f)=minstdarray;
%     bar3(minstdarray)
    % for m=1:tlen
    %     figure
    %     histogram(minstdarray(m,:),40)
    % end
    %         switch min([phi1_std,phi2_std,phi1_alt_std,phi2_alt_std])
    %             case phi1_std
    %                 calcR_all=calcR1_all;
    %             case phi2_std
    %                 calcR_all=calcR2_all;
    %             case phi1_alt_std
    %                 calcR_all=calcR3_all;
    %             case phi2_alt_std
    %                 calcR_all=calcR4_all;
    %         end
    %         axs=zeros(3,ptsnum);
    %         angs=zeros(1,ptsnum);
    %         for n=1:ptsnum
    %             axang=rotm2axang(calcR_all{n});
    %             axs(:,n)=axang(1:3)';
    %             angs(n)=axang(4);
    %         end
    %         ax=mean(axs,2);
    %         ax=ax/norm(ax);
    %         ang=mean(angs);
    %         ax_cross=[0,-ax(3),ax(2);ax(3),0,-ax(1);-ax(2),ax(1),0];
    %         ax_cross2=ax_cross*ax_cross;
    %         calcR=eye(3)+sin(ang)*ax_cross+(1-cos(ang))*ax_cross2;
    %         error=acos((trace(Rinv*calcR)-1)/2);
    %         errormap(a,b,c,d)=error;
    % for m=1:tlen
    %     figure
    %     histogram(minstdarray(m,:),40)
    % end
end

for m=1:tlen
    figure
    bar3(minstdmultiarray(:,:,m))
end
% for m=1:tlen
%     figure
%     histogram(minstdmultiarray(m,:,34),50)
% end
stdsumlist=-ones(1,tlen);
for m=1:tlen
    figure
    for n=1:tlen
        stdsumlist(n)=sum(minstdmultiarray(n,:,m));
    end
    plot(linspace(0,tlen-1,tlen),stdsumlist)
    xlabel('epipole error')
    ylabel('sum of std')
end