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
e12=[e1Theta,e1Phi,e2Theta,e2Phi].*180./pi;
epi1_sph=[e1Theta;e1Phi];
epi2_sph=[e2Theta;e2Phi];
tol=1e-2;
maxiter=10000;
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
counter_GDRE=1;
e1counter=1;
e2counter=1;
e1Theta=acos(1/norm(e1));
e2Theta=acos(1/norm(e2));
e1cosTheta=cos(e1Theta);
e1sinTheta=sin(e1Theta);
e2cosTheta=cos(e2Theta);
e2sinTheta=sin(e2Theta);
e1Phi=atan(e1(2)/e1(1));
% randr=normrnd(0,0.008,[2,ptsnum]); % Introduce noise.
% randtheta=rand(2,ptsnum).*2.*pi;
% p1noise=[randr(1,:)+cos(randtheta(1,:));randr(1,:)+sin(randtheta(1,:));zeros(1,ptsnum)];
% p2noise=[randr(2,:)+cos(randtheta(2,:));randr(2,:)+sin(randtheta(2,:));zeros(1,ptsnum)];
% p1(:,1:ptsnum)=p1(:,1:ptsnum)+p1noise;
% p2(:,1:ptsnum)=p2(:,1:ptsnum)+p2noise;
fun_GDRE=@(x)epipole_corrs_to_RE(x,p1,p2,ptsnum,'m');
lb=[0,-pi,0,-pi];
ub=[pi/2-0.001,pi,pi/2-0.001,pi];
options=optimset('MaxFunEvals',20000,'MaxIter',100000,'TolFun',1e-4,'TolX',1e-5);
total_cell={};
if e1(1)<0
    e1Phi=e1Phi+pi;
end
e2Phi=atan(e2(2)/e2(1));
if e2(1)<0
    e2Phi=e2Phi+pi;
end
epi1_sph=[e1Theta;e1Phi];
% e1s_sph(:,1)=epi1_sph;
% e2s_sph(:,1)=epi2_sph;
e12s_sph_ordered=zeros(size(e1s_sph,2)*size(e2s_sph,2),5);
e12s_sph_counter=1;
top_frac=0.1; % select top_frac of pairs that give the least RE.
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
    e1_sph=e1s_sph(:,m);
    for n=1:size(e2s_sph,2)
        e2_sph=e2s_sph(:,n);
        x=[e1_sph',e2_sph'];
        min_val_initial=norm(fun_GDRE(x)).*K(1);
        e12s_sph_ordered(e12s_sph_counter,1:5)=[x,min_val_initial];
        e12s_sph_counter=e12s_sph_counter+1;
    end
end
e12s_sph_ordered=sortrows(e12s_sph_ordered,5);
e12s_sph_ordered=e12s_sph_ordered(1:round(size(e12s_sph_ordered,1)*top_frac),:);
for n=1:size(e12s_sph_ordered,1)
    e1_sph=e12s_sph_ordered(n,1:2);
    e2_sph=e12s_sph_ordered(n,3:4);
    x=[e1_sph,e2_sph];
    min_val_initial=fun_GDRE(x);
    [xopt,fopt]=lsqnonlin(fun_GDRE,x,lb,ub);
    total_cell{counter_GDRE,1}=[e1_sph,e2_sph].*180./pi;
    total_cell{counter_GDRE,2}=xopt.*180./pi;
    total_cell{counter_GDRE,3}=min_val_initial;
    total_cell{counter_GDRE,4}=fopt;
    total_cell{counter_GDRE,5}=norm(fopt).*K(1);
    counter_GDRE=counter_GDRE+1;
end
total_cell=sortrows(total_cell,5);
total_array=reshape(cell2mat(total_cell(1:12,2)'),[4,12]);
total_array=[e12',total_array];