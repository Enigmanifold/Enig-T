rng(10)
run('generate_RT.m');
e1Theta=acos(1/norm(e1));
e2Theta=acos(1/norm(e2));
e1cosTheta=cos(e1Theta);
e1sinTheta=sin(e1Theta);
e2cosTheta=cos(e2Theta);
e2sinTheta=sin(e2Theta);
e1Phi=atan(e1(2)/e1(1));
if e1(1)<0
    e1Phi=e1Phi+pi;
    if e1Phi>pi
        e1Phi=e1Phi-2*pi;
    end
end
e2Phi=atan(e2(2)/e2(1));
if e2(1)<0
    e2Phi=e2Phi+pi;
    if e2Phi>pi
        e2Phi=e2Phi-2*pi;
    end
end
epi1_sph=[e1Theta,e1Phi]';
epi2_sph=[e2Theta,e2Phi]';
ptsnum=30;
maxfunevals=10000;
maxiter=20000;
densification2=3;
% randr_sigmas=[0,0.001,0.002];
randr_sigmas=[0,0.001,0.002,0.003,0.004,0.005,0.006,0.008,0.01];
trials=length(randr_sigmas);
psruns=50;
GD_xs=-ones(trials+1,4);
PS_xs=-ones(trials+1,4);
GD_xs(1,2:4)=[epi2_sph',0].*180./pi;
PS_xs(1,2:4)=[epi2_sph',0].*180./pi;
sorted_xstds={};
counter1=1;
counter2=1;
counter3=1;
noisy_p1s={};
noisy_p2s={};
total_cells={};
for trial=1:trials
    randr=normrnd(0,randr_sigmas(trial),[2,ptsnum]); % Introduce noise.
    randtheta=rand(2,ptsnum).*2.*pi;
    p1noise=[randr(1,:).*cos(randtheta(1,:));randr(1,:).*sin(randtheta(1,:));zeros(1,ptsnum)];
    p2noise=[randr(2,:).*cos(randtheta(2,:));randr(2,:).*sin(randtheta(2,:));zeros(1,ptsnum)];
    noisy_p1(:,1:ptsnum)=p1(:,1:ptsnum)+p1noise;
    noisy_p2(:,1:ptsnum)=p2(:,1:ptsnum)+p2noise;
    noisy_p1s{counter1}=noisy_p1;
    noisy_p2s{counter1}=noisy_p2;
    counter1=counter1+1;
end
for trial=1:trials
    p1_with_noise=noisy_p1s{trial};
    p2_with_noise=noisy_p2s{trial};
    GD_xs(trial+1,1)=randr_sigmas(trial);
    [GD_x,total_cell]=GD_image_epi1known(epi1_sph,p1_with_noise,p2_with_noise,ptsnum,maxfunevals,maxiter,densification2);
    GD_xs(trial+1,2:3)=GD_x*180./pi;
    GD_xs(trial+1,4)=acos(sin(epi2_sph(1))*sin(GD_x(1))+cos(epi2_sph(1))*cos(GD_x(1))*cos(epi2_sph(2)-GD_x(2))).*180./pi;
    total_cells{counter2}=total_cell;
    counter2=counter2+1;
end    
for trial=1:trials
    p1_with_noise=noisy_p1s{trial};
    p2_with_noise=noisy_p2s{trial};
    PS_xs(trial+1,1)=randr_sigmas(trial);
    [PS_x,sorted_xstd]=PS_image_epi1known(epi1_sph,p1_with_noise,p2_with_noise,ptsnum,psruns);
    PS_xs(trial+1,2:3)=PS_x.*180./pi;
    PS_xs(trial+1,4)=acos(sin(epi2_sph(1))*sin(PS_x(1))+cos(epi2_sph(1))*cos(PS_x(1))*cos(epi2_sph(2)-PS_x(2))).*180./pi;
    sorted_xstds{counter3}=sorted_xstd;
    counter3=counter3+1;
end