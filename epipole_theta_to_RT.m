%% Getting epipoles and theta.
run('generate_RT.m');
e1=R'*T;
e1=e1/e1(3);
e2=T/T(3);
ptsnum=size(DPoints,2);
ptsnum=2000;
p1_to_e1=p1(:,1:ptsnum)-repmat(e1,1,ptsnum);
p2_to_e2=p2(:,1:ptsnum)-repmat(e2,1,ptsnum);
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
% v2n=zeros(3,ptsnum);
un=zeros(1,ptsnum);
uvn=zeros(1,ptsnum);
% a=zeros(1,ptsnum);
% b=zeros(1,ptsnum);
% c=zeros(1,ptsnum);
vn_alt=zeros(3,ptsnum);
% v2n_alt=zeros(3,ptsnum);
un_alt=zeros(1,ptsnum);
uvn_alt=zeros(1,ptsnum);
% a_alt=zeros(1,ptsnum);
% b_alt=zeros(1,ptsnum);
% c_alt=zeros(1,ptsnum);
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
% for m=1:ptsnum
%     a(m)=-u0_bar(:,m)'*v2n(:,m);
%     b(m)=u0_bar(:,m)'*vn(:,m);
%     c(m)=-un(m)-u0_bar(:,m)'*v2n(:,m);
%     a_alt(m)=-u0_bar(:,m)'*v2n_alt(:,m);
%     b_alt(m)=u0_bar(:,m)'*vn_alt(:,m);
%     c_alt(m)=-un(m)-u0_bar(:,m)'*v2n_alt(:,m);
% end
% abnorm=sqrt(a.*a+b.*b);
% abnorm_alt=sqrt(a_alt.*a_alt+b_alt.*b_alt);
% cosgamma=a./abnorm;
% singamma=b./abnorm;
% cosgamma_alt=a_alt./abnorm_alt;
% singamma_alt=b_alt./abnorm_alt;
% gamma=acos(cosgamma);
% gamma_alt=acos(cosgamma_alt);
% for m=1:ptsnum
%     if singamma(m)<0
%         gamma(m)=2*pi-gamma(m);
%     end
%     if singamma_alt(m)<0
%         gamma_alt(m)=2*pi-gamma_alt(m);
%     end
% end
% phi1=gamma+repmat(0.5*pi,1,ptsnum);
% phi2=gamma+repmat(1.5*pi,1,ptsnum);
% phi1_alt=gamma_alt+repmat(0.5*pi,1,ptsnum);
% phi2_alt=gamma_alt+repmat(1.5*pi,1,ptsnum);
% for m=1:ptsnum
%     if phi2(m)>2*pi
%         phi2(m)=phi2(m)-2*pi;
%     end
%     if phi2_alt(m)>2*pi
%         phi2_alt(m)=phi2_alt(m)-2*pi;
%     end
% end
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
calcR1=cell(1,ptsnum);
calcR2=cell(1,ptsnum);
calcR3=cell(1,ptsnum);
calcR4=cell(1,ptsnum);
for m=1:ptsnum
    calcR1{m}=R2_1{m}*R1;
    calcR2{m}=R2_2{m}*R1;
    calcR3{m}=R2_3{m}*R1_alt;
    calcR4{m}=R2_4{m}*R1_alt;
end

error=-ones(1,ptsnum);
ind=-ones(1,ptsnum);
Rinv=R';
for m=1:ptsnum
    if isreal(calcR1{m})
        R1error=acos((trace(Rinv*calcR1{m})-1)/2);
    else
        R1error=99;
    end
    if isreal(calcR2{m})
        R2error=acos((trace(Rinv*calcR2{m})-1)/2);
    else
        R2error=99;
    end
    if isreal(calcR3{m})
        R3error=acos((trace(Rinv*calcR3{m})-1)/2);
    else
        R3error=99;
    end
    if isreal(calcR4{m})
        R4error=acos((trace(Rinv*calcR4{m})-1)/2);
    else
        R4error=99;
    end
    [error(m),ind(m)]=min([R1error,R2error,R3error,R4error]);
end
if min(ind)==max(ind)
    disp(strcat('From calcR',num2str(min(ind))));
    switch min(ind)
        case 1
            calcR=calcR1{1};
        case 2
            calcR=calcR2{1};
        case 3
            calcR=calcR3{1};
        case 4
            calcR=calcR4{1};
    end
else
    disp('Indeterminate R!');
end
% rotate_calcR=pi_around_T*calcR;
% epi1=rotate_calcR'*T;
% epi1=epi1/epi1(3);

check_depth1=zeros(1,ptsnum);
check_depth2=zeros(1,ptsnum);
check_depth3=zeros(1,ptsnum);
check_depth4=zeros(1,ptsnum);
for m=1:ptsnum
    check_depth1(m)=v2'*cross(u0_bar(:,m),cross(calcR1{m}*u0(:,m),v2));
    check_depth2(m)=v2'*cross(u0_bar(:,m),cross(calcR2{m}*u0(:,m),v2));
    check_depth3(m)=v2'*cross(u0_bar(:,m),cross(calcR3{m}*u0(:,m),v2));
    check_depth4(m)=v2'*cross(u0_bar(:,m),cross(calcR4{m}*u0(:,m),v2));
end
disp(strcat('Total number of points: ',num2str(ptsnum)));
disp(strcat('Number of positive depths for calcR1 is ',num2str(sum(check_depth1<0))));
disp(strcat('Number of positive depths for calcR2 is ',num2str(sum(check_depth2<0))));
disp(strcat('Number of positive depths for calcR3 is ',num2str(sum(check_depth3<0))));
disp(strcat('Number of positive depths for calcR4 is ',num2str(sum(check_depth4<0))));
% for m=1:ptsnum
%     vn(:,m)=v1_cross*N2_array(:,m);
%     vn_alt(:,m)=v1_cross*N2_array_alt(:,m);
%     v2n(:,m)=v1_cross2*N2_array(:,m);
%     v2n_alt(:,m)=v1_cross2*N2_array_alt(:,m);
%     un(:,m)=u0_bar(:,m)'*N2_array(:,m);
%     un_alt(:,m)=u0_bar(:,m)'*N2_array_alt(:,m);
% end
% for m=1:ptsnum
%     a(m)=-u0_bar(:,m)'*v2n(:,m);
%     b(m)=u0_bar(:,m)'*vn(:,m);
%     c(m)=-un(m)-u0_bar(:,m)'*v2n(:,m);
%     a_alt(m)=-u0_bar(:,m)'*v2n_alt(:,m);
%     b_alt(m)=u0_bar(:,m)'*vn_alt(:,m);
%     c_alt(m)=-un(m)-u0_bar(:,m)'*v2n_alt(:,m);
% end
% abnorm=sqrt(a.*a+b.*b);
% abnorm_alt=sqrt(a_alt.*a_alt+b_alt.*b_alt);
% cosgamma=a./abnorm;
% singamma=b./abnorm;
% cosgamma_alt=a_alt./abnorm_alt;
% singamma_alt=b_alt./abnorm_alt;
% gamma=acos(cosgamma);
% gamma_alt=acos(cosgamma_alt);
% for m=1:ptsnum
%     if singamma(m)<0
%         gamma(m)=2*pi-gamma(m);
%     end
%     if singamma_alt(m)<0
%         gamma_alt(m)=2*pi-gamma_alt(m);
%     end
% end
% phi1=gamma+acos(c./abnorm);
% phi2=gamma+repmat(pi,1,ptsnum)-acos(c./abnorm);
% phi1_alt=gamma_alt+acos(c_alt./abnorm_alt);
% phi2_alt=gamma_alt+repmat(pi,1,ptsnum)-acos(c_alt./abnorm_alt);
% sinphi1=sin(phi1);
% cosphi1=cos(phi1);
% sinphi2=sin(phi2);
% cosphi2=cos(phi2);
% sinphi1_alt=sin(phi1_alt);
% cosphi1_alt=cos(phi1_alt);
% sinphi2_alt=sin(phi2_alt);
% cosphi2_alt=cos(phi2_alt);
% R2_1=cell(1,ptsnum);
% R2_2=cell(1,ptsnum);
% R2_3=cell(1,ptsnum);
% R2_4=cell(1,ptsnum);
% for m=1:ptsnum
%     R2_1{m}=eye(3)+sinphi1(m)*v1_cross+(1-cosphi1(m))*v1_cross2;
%     R2_2{m}=eye(3)+sinphi2(m)*v1_cross+(1-cosphi2(m))*v1_cross2;
%     R2_3{m}=eye(3)+sinphi1_alt(m)*v1_cross+(1-cosphi1_alt(m))*v1_cross2;
%     R2_4{m}=eye(3)+sinphi2_alt(m)*v1_cross+(1-cosphi2_alt(m))*v1_cross2;
% end
% calcR1=cell(1,ptsnum);
% calcR2=cell(1,ptsnum);
% calcR3=cell(1,ptsnum);
% calcR4=cell(1,ptsnum);
% for m=1:ptsnum
%     calcR1{m}=R2_1{m}*R1;
%     calcR2{m}=R2_2{m}*R1;
%     calcR3{m}=R2_3{m}*R1;
%     calcR4{m}=R2_4{m}*R1;
% end
% theta_diff=theta2-theta1;
% for m=1:ptsnum
%     if theta_diff(m)>0
%         theta_diff(m)=theta_diff(m)-2*pi;
%     end
% end
% plot(linspace(0,ptsnum-1,ptsnum),theta_diff)
% rerror=acos((trace(inv(calcR1{1})*calcR2)-1)/2);