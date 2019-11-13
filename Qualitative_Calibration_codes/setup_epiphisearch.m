% run('generate_RT.m');

e1Theta=acos(1/norm(e1));
e2Theta=acos(1/norm(e2));
e1cosTheta=cos(e1Theta);
e1sinTheta=sin(e1Theta);
e2cosTheta=cos(e2Theta);
e2sinTheta=sin(e2Theta);
e1Phi=atan(e1(2)/e1(1));
phi_div=10;
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
e12=[e1Theta,e1Phi,e2Theta,e2Phi].*180./pi;
epi1_sph=[e1Theta;e1Phi];
epi2_sph=[e2Theta;e2Phi];
epi12_sphgt=[epi1_sph',epi2_sph'];
ptsnum=30;
e1s_carte=IcosahedronMesh;
e2s_carte=IcosahedronMesh;
k1=2;
k2=2;
e1s_carte=SubdivideSphericalMesh(e1s_carte,k1);
e2s_carte=SubdivideSphericalMesh(e2s_carte,k2);
e1s_sph=[];
e2s_sph=[];
e1counter=1;
e2counter=1;
e1Theta=acos(1/norm(e1));
e2Theta=acos(1/norm(e2));
e1cosTheta=cos(e1Theta);
e1sinTheta=sin(e1Theta);
e2cosTheta=cos(e2Theta);
e2sinTheta=sin(e2Theta);
e1Phi=atan(e1(2)/e1(1));
% randr=normrnd(0,0.01,[2,ptsnum]); % Introduce noise.
% randtheta=rand(2,ptsnum).*2.*pi;
% p1noise=[randr(1,:).*cos(randtheta(1,:));randr(1,:).*sin(randtheta(1,:));zeros(1,ptsnum)];
% p2noise=[randr(2,:).*cos(randtheta(2,:));randr(2,:).*sin(randtheta(2,:));zeros(1,ptsnum)];
% p1(:,1:ptsnum)=p1(:,1:ptsnum)+p1noise;
% p2(:,1:ptsnum)=p2(:,1:ptsnum)+p2noise;
if e1(1)<0
    e1Phi=e1Phi+pi;
end
e2Phi=atan(e2(2)/e2(1));
if e2(1)<0
    e2Phi=e2Phi+pi;
end
phi_list=linspace(0,2*pi*(1-1/phi_div),phi_div);
for m=1:size(e1s_carte.Points,1)
    if e1s_carte.Points(m,3)>0
        if e1s_carte.Points(m,1)~=0
            sph_coords=[acos(e1s_carte.Points(m,3));atan(e1s_carte.Points(m,2)/e1s_carte.Points(m,1))];
        else
            sph_coords=[acos(e1s_carte.Points(m,3));pi/2];
        end
        if e1s_carte.Points(m,1)<0 && e1s_carte.Points(m,2)>=0
            sph_coords(2)=sph_coords(2)+pi;
        end
        if e1s_carte.Points(m,1)<=0 && e1s_carte.Points(m,2)<0
            sph_coords(2)=sph_coords(2)-pi;
        end
    e1counter=e1counter+1;
    end
    e1s_sph(:,e1counter)=sph_coords;
end
for m=1:size(e2s_carte.Points,1)
    if e2s_carte.Points(m,3)>0
        if e2s_carte.Points(m,1)~=0
            sph_coords=[acos(e2s_carte.Points(m,3));atan(e2s_carte.Points(m,2)/e2s_carte.Points(m,1))];
        else
            sph_coords=[acos(e2s_carte.Points(m,3));pi/2];
        end
        if e2s_carte.Points(m,1)<0 && e2s_carte.Points(m,2)>=0
            sph_coords(2)=sph_coords(2)+pi;
        end
        if e2s_carte.Points(m,1)<=0 && e2s_carte.Points(m,2)<0
            sph_coords(2)=sph_coords(2)-pi;
        end
        e2counter=e2counter+1;
    end
    e2s_sph(:,e2counter)=sph_coords;
end
e1s_sphsize=size(e1s_sph,2);
e2s_sphsize=size(e2s_sph,2);
epiphi=-ones(e1s_sphsize*e2s_sphsize*phi_div,5);
epiphi_counter=1;
for m=1:e1s_sphsize
    for n=1:e2s_sphsize
        for o=1:phi_div
            epiphi(epiphi_counter,:)=[e1s_sph(:,m)',e2s_sph(:,n)',phi_list(o)];
            epiphi_counter=epiphi_counter+1;
        end
    end
end