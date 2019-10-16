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
e12=[e1Theta,e1Phi,e2Theta,e2Phi];
ptsnum=30;
randr_sigmas=[0,0.001,0.002,0.003,0.004,0.005,0.007,0.01];
trials=length(randr_sigmas);
psruns=200;
final_xs=-ones(trials+1,8);
final_xs(1,2:end)=[e12,0,0,0];
sorted_xstds={};
counter1=1;
counter2=1;
noisy_p1s={};
noisy_p2s={};
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
    final_xs(trial+1,1)=randr_sigmas(trial);
    [final_x,sorted_xstd]=PS_image(p1_with_noise,p2_with_noise,ptsnum,psruns);
    final_xs(trial+1,2:5)=final_x;
    final_xs(trial+1,6)=acos(sin(e12(1))*sin(final_x(1))+cos(e12(1))*cos(final_x(1))*cos(e12(2)-final_x(2))).*180./pi;
    final_xs(trial+1,7)=acos(sin(e12(3))*sin(final_x(3))+cos(e12(3))*cos(final_x(3))*cos(e12(4)-final_x(4))).*180./pi;
    final_xs(trial+1,8)=final_xs(trial+1,6)+final_xs(trial+1,7);
    sorted_xstds{counter2}=sorted_xstd;
    counter2=counter2+1;
end