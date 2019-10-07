%% Errors due to increasing constant distances.
run('generate_RT.m');
Rinv=R';
e1=R'*T;
e1=e1/e1(3);
e2=T/T(3);
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
t=linspace(-50,0,30).*0.1;
errs=[flip(-exp(t)),0,exp(t)];
tlen=length(errs);
samplen=100;
minstdarray=-ones(tlen,samplen);
minstdmultiarray=-ones(tlen,samplen,tlen);
total_cell={};
counter=1;
theta_overflow=0;
for f=1:tlen
    rndangles1=errs(f).*rand.*2.*pi;
    rndepi1err_sph=[errs(f);rndangles1];
    epi1_samples_sph=e1_sph+rndepi1err_sph;
    for a=1:tlen
        err=errs(a);
%         rndangles1=rand(1,samplen).*2.*pi;
        rndangles2=linspace(-err,err,samplen).*pi;
%         rndepi1err=[errs(f).*cos(rndangles1);errs(f).*sin(rndangles1);zeros(1,samplen)];
        rndepi2err_sph=[ones(1,samplen).*err;rndangles2;];
        epi2_samples_sph=repmat(e2_sph,1,samplen)+rndepi2err_sph;
        for b=1:samplen
            epi1_sph=epi1_samples_sph;
            epi2_sph=epi2_samples_sph(:,b);
            if epi1_sph(1)>=pi/2 || epi2_sph(1)>=pi/2
                theta_overflow=theta_overflow+1;
                continue
            end
            epi1=[sin(epi1_sph(1))*cos(epi1_sph(2))*sec(epi1_sph(1));sin(epi1_sph(1))*sin(epi1_sph(2))*sec(epi1_sph(1));1];
            epi2=[sin(epi2_sph(1))*cos(epi2_sph(2))*sec(epi2_sph(1));sin(epi2_sph(1))*sin(epi2_sph(2))*sec(epi2_sph(1));1];
            [phi1,phi2,phi1_alt,phi2_alt]=calc_phi(p1(:,1:ptsnum),p2(:,1:ptsnum),epi1,epi2,ptsnum);
            phis=[phi1;phi2;phi1_alt;phi2_alt];
            phis_c0=-ones(4,ptsnum);
            for c=1:4               
                for d=1:ptsnum
                    [~,phi2_c0_ind]=min([abs(phis(c,d)),abs(phis(c,d)-pi),abs(phis(c,d)-2*pi)]);
                    switch phi2_c0_ind
                        case 1
                            phis_c0(c,d)=phis(c,d);
                        case 2
                            phis_c0(c,d)=phis(c,d)-pi;
                        case 3
                            phis_c0(c,d)=phis(c,d)-2*pi;
                    end
                end
            end
            phi1_c0std=std(phis_c0(1,:));
            phi2_c0std=std(phis_c0(2,:));
            phi3_c0std=std(phis_c0(3,:));
            phi4_c0std=std(phis_c0(4,:));
            [minc0std,minc0stdind]=min([phi1_c0std,phi2_c0std,phi3_c0std,phi4_c0std]);
            chosen_phis=phis(minc0stdind,:);
            minstdarray(a,b)=minc0std;
        end
    end
    minstdmultiarray(:,:,f)=minstdarray;
end
for m=1:tlen
    figure
    bar3(minstdmultiarray(:,:,m))
end       
% stdsumlist=-ones(1,tlen);
% for m=1:tlen
%     figure
%     for n=1:tlen
%         stdsumlist(n)=sum(minstdmultiarray(n,:,m));
%     end
%     plot(linspace(0,tlen-1,tlen),stdsumlist)
%     xlabel('epipole error')
%     ylabel('sum of std')
% end    