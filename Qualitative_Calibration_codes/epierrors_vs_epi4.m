%% Errors due to increasing constant distances.
run('generate_RT.m');
Rinv=R';
e1=R'*T;
e1=e1/e1(3);
e2=T/T(3);
ptsnum=size(DPoints,2);
ptsnum=30;
t=linspace(-50,50,100).*0.1;
tlen=length(t);
errs=exp(t);
errs(1)=0;
samplen=100;
minstdarray=-ones(tlen,samplen);
minstdmultiarray=-ones(tlen,samplen,tlen);
total_cell={};
counter=1;
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
            epi1=epi1_samples(:,b);
            epi2=epi2_samples(:,b);
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