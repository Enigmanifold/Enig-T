run('generate_RT.m');
tol=1e-2;
maxiter=100000;
angle_div=4;
alpha=0.001;
beta = linspace(0,2*pi*(1-1/angle_div),angle_div);
ptsnum=30;
diff_amplifier=1000000;
e1s_carte=IcosahedronMesh;
e2s_carte=IcosahedronMesh;
k1=2;
k2=2;
e1s_carte=SubdivideSphericalMesh(e1s_carte,k1);
e2s_carte=SubdivideSphericalMesh(e2s_carte,k2);
e1s_sph=[];
e2s_sph=[];
counter=1;
e1counter=1;
e2counter=1;
e1Theta=acos(1/norm(e1));
e2Theta=acos(1/norm(e2));
e1cosTheta=cos(e1Theta);
e1sinTheta=sin(e1Theta);
e2cosTheta=cos(e2Theta);
e2sinTheta=sin(e2Theta);
e1Phi=atan(e1(2)/e1(1));
% randr=normrnd(0,0.004,[2,ptsnum]); % Introduce noise.
% randtheta=rand(2,ptsnum).*2.*pi;
% p1noise=[randr(1,:).*cos(randtheta(1,:));randr(1,:).*sin(randtheta(1,:));zeros(1,ptsnum)];
% p2noise=[randr(2,:).*cos(randtheta(2,:));randr(2,:).*sin(randtheta(2,:));zeros(1,ptsnum)];
% p1(:,1:ptsnum)=p1(:,1:ptsnum)+p1noise;
% p2(:,1:ptsnum)=p2(:,1:ptsnum)+p2noise;
fun=@(x)calc_phi3(x,p1,p2,ptsnum,diff_amplifier);
options=optimset('MaxFunEvals',10000,'MaxIter',maxiter,'TolFun',1e-3);
total_cell={};
if e1(1)<0
    e1Phi=e1Phi+pi;
end
e2Phi=atan(e2(2)/e2(1));
if e2(1)<0
    e2Phi=e2Phi+pi;
end
e1_sph=[e1Theta;e1Phi];
e2_sph=[e2Theta;e2Phi];
e1s_sph(:,1)=e1_sph;
e2s_sph(:,1)=e2_sph;
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
for m=1:size(e1s_sph,2)
    for n=1:size(e2s_sph,2)
        epi1_sph=e1s_sph(:,m);
        epi2_sph=e2s_sph(:,n);
        x=[epi1_sph;epi2_sph]';
        min_val_initial=fun(x);
        [xopt,fopt]=fminsearch(fun,x,options);
        total_cell{counter,1}=epi1_sph';
        total_cell{counter,2}=epi2_sph';
        total_cell{counter,3}=xopt(1:2);
        total_cell{counter,4}=xopt(3:4);
        total_cell{counter,5}=min_val_initial;
        total_cell{counter,6}=fopt;
        counter=counter+1;
    end
end
total_cell=sortrows(total_cell,6);
total_array=reshape(cell2mat(total_cell(1:12,4)'),[2,12]).*180./pi;